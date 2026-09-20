"""跨字段语义校验：位点越界、重复位点、酶名称重复等。

校验失败统一抛出 :class:`DigestValidationError`，由主应用转换为
统一结构的 HTTP 422 响应（与 Pydantic 字段校验同一响应格式）。
"""
from fastapi import HTTPException

from app.config import settings
from app.models import DigestRequest, IdentifyRequest, PanelRequest


class DigestValidationError(HTTPException):
    def __init__(self, errors: list[dict]):
        super().__init__(status_code=422, detail={"code": "validation_failed", "errors": errors})


def _err(loc: list[str], message: str, code: str = "invalid") -> dict:
    return {"loc": loc, "msg": message, "type": code}


def _check_request_size(errors: list[dict], req: DigestRequest | PanelRequest,
                        length: int) -> None:
    """最小可见片段与酶数量上限的公共校验。"""
    if req.min_visible_fragment >= length:
        errors.append(_err(
            ["min_visible_fragment"],
            f"最小可见片段长度（{req.min_visible_fragment}）必须小于质粒长度（{length}）",
            "min_visible_too_large",
        ))

    if len(req.enzymes) > settings.MAX_ENZYMES:
        errors.append(_err(
            ["enzymes"],
            f"酶数量（{len(req.enzymes)}）超过允许上限 {settings.MAX_ENZYMES}",
            "too_many_enzymes",
        ))


def _check_sites(errors: list[dict], loc: list, label: str, sites: list[int],
                 length: int) -> None:
    """校验一组切点的越界与重复；label 为中文场景描述（如"酶 EcoRI 在参考型中"）。"""
    out_of_range = [s for s in sites if s < 0 or s >= length]
    if out_of_range:
        errors.append(_err(
            loc,
            f"{label}存在越界切点 {out_of_range}：合法范围为 0 ~ {length - 1}",
            "site_out_of_range",
        ))

    unique_sites = set(sites)
    if len(unique_sites) != len(sites):
        dupes = sorted({s for s in sites if sites.count(s) > 1})
        errors.append(_err(loc, f"{label}存在重复切点 {dupes}", "duplicate_site"))


def validate_request(req: DigestRequest) -> DigestRequest:
    errors: list[dict] = []
    length = req.plasmid_length
    _check_request_size(errors, req, length)

    seen_names: set[str] = set()
    for i, enzyme in enumerate(req.enzymes):
        name = enzyme.name
        if name in seen_names:
            errors.append(_err(
                ["enzymes", i, "name"],
                f"酶名称重复：{name}",
                "duplicate_enzyme_name",
            ))
        seen_names.add(name)

        for field in ("reference_sites", "variant_sites"):
            sites = getattr(enzyme, field)
            sample = "参考型" if field == "reference_sites" else "变异型"
            _check_sites(errors, ["enzymes", i, field], f"酶 {name} 在{sample}中",
                         sites, length)

    if errors:
        raise DigestValidationError(errors)
    return req


def validate_panel_request(req: PanelRequest) -> PanelRequest:
    """panel 接口的语义校验：位点越界 / 重复、酶名重复、空变异名、规模上限。"""
    errors: list[dict] = []
    length = req.plasmid_length
    _check_request_size(errors, req, length)

    seen_names: set[str] = set()
    for i, enzyme in enumerate(req.enzymes):
        name = enzyme.name
        if name in seen_names:
            errors.append(_err(
                ["enzymes", i, "name"],
                f"酶名称重复：{name}",
                "duplicate_enzyme_name",
            ))
        seen_names.add(name)

        _check_sites(errors, ["enzymes", i, "reference_sites"],
                     f"酶 {name} 在参考型中", enzyme.reference_sites, length)
        for variant, sites in enzyme.variant_sites.items():
            if not variant:
                errors.append(_err(
                    ["enzymes", i, "variant_sites"],
                    f"酶 {name} 的 variant_sites 中存在空变异名",
                    "empty_variant_name",
                ))
            _check_sites(errors, ["enzymes", i, "variant_sites", variant],
                         f"酶 {name} 在变异 {variant} 中", sites, length)

    if errors:
        raise DigestValidationError(errors)
    return req


def validate_identify_request(req: IdentifyRequest) -> IdentifyRequest:
    """identify 接口的语义校验：条带列表规模上限。

    观测为空、条带为负、容差非法等基础字段约束由 Pydantic 模型直接拦截，
    与此处一样输出统一结构的 422 响应。
    """
    errors: list[dict] = []
    for field in ("reference_bands", "variant_bands", "observed_bands"):
        bands = getattr(req, field)
        if len(bands) > settings.MAX_BANDS:
            errors.append(_err(
                [field],
                f"条带数量（{len(bands)}）超过允许上限 {settings.MAX_BANDS}",
                "too_many_bands",
            ))
    if errors:
        raise DigestValidationError(errors)
    return req
