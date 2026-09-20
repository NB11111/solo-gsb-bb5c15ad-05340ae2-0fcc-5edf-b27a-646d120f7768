"""端到端 API 测试：最小酶对面板（POST /api/v1/digest/panel）及精确求解器。"""
from fastapi.testclient import TestClient

from app.main import app
from app.panel import PairCoverage, select_minimum_panel

client = TestClient(app)


def _enzyme(name, ref, variants=None):
    enzyme = {"name": name, "reference_sites": ref}
    if variants is not None:
        enzyme["variant_sites"] = variants
    return enzyme


def _payload(enzymes, length=1000, min_visible=50, tolerance=20):
    return {
        "plasmid_length": length,
        "min_visible_fragment": min_visible,
        "band_resolution_tolerance": tolerance,
        "enzymes": enzymes,
    }


def _cand(key, mask, bands):
    return PairCoverage(key.split("+"), key, mask, bands, {})


# ---------------------------------------------------------------- 最小面板

def test_panel_minimum_with_band_tie_break():
    # L=1000、最小可见=50、容差=20，手工推演：
    #   A+B：仅覆盖 v1（2 条可区分条带）；A+C：仅覆盖 v2（2 条）；
    #   A+D：仅覆盖 v1（4 条）；B+C / B+D / C+D：参考型仅 1 个切点，不可用。
    # 最小面板规模为 2：{A+B,A+C} 共 4 条 < {A+D,A+C} 共 6 条 -> 选后者。
    payload = _payload([
        _enzyme("A", [100, 600], {"v1": [100, 600], "v2": [100, 600]}),
        _enzyme("B", [300], {"v1": [300, 850], "v2": [300]}),
        _enzyme("C", [300], {"v1": [300], "v2": [300, 850]}),
        _enzyme("D", [300], {"v1": [300, 850, 900], "v2": [300]}),
    ])
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 200
    body = res.json()

    assert body["variants"] == ["v1", "v2"]
    assert body["total_enzyme_pairs"] == 6
    assert body["uncovered_variants"] == []

    pairs = body["selected_pairs"]
    assert [p["pair_key"] for p in pairs] == ["A+C", "A+D"]
    assert pairs[0]["covers_variants"] == ["v2"]
    assert pairs[0]["distinguishable_band_count"] == 2
    assert pairs[1]["covers_variants"] == ["v1"]
    assert pairs[1]["distinguishable_band_count"] == 4

    coverage = body["variant_coverage"]
    assert [c["variant"] for c in coverage] == ["v1", "v2"]
    assert coverage[0]["pair_key"] == "A+D"
    assert coverage[0]["distinguishable_band_count"] == 4
    assert len(coverage[0]["differences"]) == 4
    assert coverage[1]["pair_key"] == "A+C"
    assert coverage[1]["distinguishable_band_count"] == 2
    # 差异依据与推荐接口同结构，且每条确实可区分
    for entry in coverage:
        for diff in entry["differences"]:
            assert diff["reference_count"] != diff["variant_count"]
            assert diff["difference_type"] in {
                "reference_only", "variant_only", "count_mismatch",
            }
            assert diff["description"]


def test_panel_single_pair_covers_all_variants():
    payload = _payload([
        _enzyme("E1", [100, 600], {"v1": [100, 600, 850], "v2": [100, 600, 850]}),
        _enzyme("E2", [300], {"v1": [300], "v2": [300]}),
    ])
    body = client.post("/api/v1/digest/panel", json=payload).json()
    assert body["uncovered_variants"] == []
    assert [p["pair_key"] for p in body["selected_pairs"]] == ["E1+E2"]
    assert body["selected_pairs"][0]["covers_variants"] == ["v1", "v2"]
    assert body["selected_pairs"][0]["distinguishable_band_count"] == 4
    assert [c["variant"] for c in body["variant_coverage"]] == ["v1", "v2"]
    assert all(c["pair_key"] == "E1+E2" for c in body["variant_coverage"])


def test_panel_uncovered_variant_returns_no_partial_panel():
    # v2 上两种酶的切点与参考型完全一致 -> 任何酶对都无法区分 v2
    payload = _payload([
        _enzyme("E1", [100, 600], {"v1": [100, 600, 850], "v2": [100, 600]}),
        _enzyme("E2", [300], {"v1": [300], "v2": [300]}),
    ])
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 200
    body = res.json()
    assert body["uncovered_variants"] == ["v2"]
    assert body["selected_pairs"] == []
    assert body["variant_coverage"] == []
    assert body["variants"] == ["v1", "v2"]


def test_panel_enzyme_without_variant_entry_treated_as_no_sites():
    # E1 未提供 v2 的切点 -> 按空列表处理；v2 上仅 E2 的 2 个切点仍可成带
    payload = _payload([
        _enzyme("E1", [100, 600], {"v1": [100, 600, 850]}),
        _enzyme("E2", [300], {"v1": [300], "v2": [300, 700]}),
    ])
    body = client.post("/api/v1/digest/panel", json=payload).json()
    assert body["variants"] == ["v1", "v2"]
    assert body["uncovered_variants"] == []
    assert [p["pair_key"] for p in body["selected_pairs"]] == ["E1+E2"]
    assert body["selected_pairs"][0]["covers_variants"] == ["v1", "v2"]
    # v1 贡献 2 条、v2 贡献 5 条可区分条带
    assert body["selected_pairs"][0]["distinguishable_band_count"] == 7


def test_panel_insufficient_cuts_pair_cannot_cover():
    # 参考型仅 1 个合并切点 -> 该酶对对 v1 不可用，v1 无其他酶对 -> 未覆盖
    payload = _payload([
        _enzyme("E1", [100], {"v1": [100, 400]}),
        _enzyme("E2", []),
    ])
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 200
    body = res.json()
    assert body["total_enzyme_pairs"] == 1
    assert body["uncovered_variants"] == ["v1"]
    assert body["selected_pairs"] == []
    assert body["variant_coverage"] == []


# ---------------------------------------------------------------- 精确求解器

def test_solver_is_exact_not_greedy():
    # 贪心（先拿覆盖最多的 G+R={0,1,3,4}）需要 3 对；最优为 E+P + F+Q 共 2 对
    candidates = [
        _cand("E+P", 0b000111, 3),
        _cand("F+Q", 0b111000, 3),
        _cand("G+R", 0b011011, 4),
    ]
    chosen = select_minimum_panel(candidates, 6)
    assert sorted(c.key for c in chosen) == ["E+P", "F+Q"]


def test_solver_tie_break_prefers_more_bands():
    # A+X 与 B+Y 覆盖集合相同，B+Y 条带更多 -> 选 B+Y
    candidates = [
        _cand("A+X", 0b01, 3),
        _cand("B+Y", 0b01, 5),
        _cand("C+Z", 0b10, 1),
    ]
    chosen = select_minimum_panel(candidates, 2)
    assert sorted(c.key for c in chosen) == ["B+Y", "C+Z"]


def test_solver_tie_break_lexicographic_pair_keys():
    # {P2,P3} 与 {P3,P4} 同为 2 对、4 条带；按排序后 pair_key 列表字典序选后者
    candidates = [
        _cand("C+X", 0b001, 1),
        _cand("D+Y", 0b110, 2),
        _cand("A+Z", 0b011, 2),
        _cand("B+W", 0b100, 2),
    ]
    chosen = select_minimum_panel(candidates, 3)
    assert sorted(c.key for c in chosen) == ["A+Z", "B+W"]


def test_solver_returns_none_when_uncoverable():
    assert select_minimum_panel([_cand("A+B", 0b01, 1)], 2) is None
    assert select_minimum_panel([], 0) == []


# ---------------------------------------------------------------- 422 校验

def _expect_panel_422(payload, code_part=None):
    res = client.post("/api/v1/digest/panel", json=payload)
    assert res.status_code == 422
    body = res.json()
    assert body["code"] == "validation_failed"
    assert body["errors"]
    if code_part:
        assert any(code_part in e["type"] for e in body["errors"])
    return body


def test_panel_422_site_out_of_range():
    payload = _payload([
        _enzyme("E1", [100, 600], {"v1": [100, 1000]}),   # 1000 越界（合法 0..999）
        _enzyme("E2", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "site_out_of_range")


def test_panel_422_negative_site_in_variant():
    payload = _payload([
        _enzyme("E1", [100], {"v1": [-1]}),
        _enzyme("E2", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "site_out_of_range")


def test_panel_422_duplicate_site():
    payload = _payload([
        _enzyme("E1", [100, 100], {"v1": [100]}),
        _enzyme("E2", [300], {"v1": [300, 300]}),
    ])
    body = _expect_panel_422(payload, "duplicate_site")
    # 参考型与变异型列表的重复切点一次性全部报告
    assert sum("duplicate_site" in e["type"] for e in body["errors"]) == 2


def test_panel_422_duplicate_enzyme_name():
    payload = _payload([
        _enzyme("E1", [100], {"v1": [100]}),
        _enzyme("E1", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "duplicate_enzyme_name")


def test_panel_422_empty_variant_name():
    payload = _payload([
        _enzyme("E1", [100], {"": [100]}),
        _enzyme("E2", [300], {"v1": [300]}),
    ])
    _expect_panel_422(payload, "empty_variant_name")


def test_panel_422_min_visible_not_smaller_than_plasmid():
    payload = _payload(
        [_enzyme("E1", [100], {"v1": [100]}), _enzyme("E2", [300], {"v1": [300]})],
        min_visible=1000,
    )
    _expect_panel_422(payload, "min_visible_too_large")


def test_panel_422_too_many_enzymes(monkeypatch):
    from app.config import settings
    monkeypatch.setattr(settings, "MAX_ENZYMES", 2)
    payload = _payload([
        _enzyme("E1", [100], {"v1": [100]}),
        _enzyme("E2", [300], {"v1": [300]}),
        _enzyme("E3", [500], {"v1": [500]}),
    ])
    _expect_panel_422(payload, "too_many_enzymes")


def test_panel_422_basic_field_constraints():
    # 质粒长度 <= 0、容差为负、酶不足 2 种、缺字段
    _expect_panel_422(_payload(
        [_enzyme("E1", [100], {"v1": [100]}), _enzyme("E2", [300], {"v1": [300]})],
        length=0,
    ))
    _expect_panel_422(_payload(
        [_enzyme("E1", [100], {"v1": [100]}), _enzyme("E2", [300], {"v1": [300]})],
        tolerance=-1,
    ))
    _expect_panel_422(_payload([_enzyme("E1", [100], {"v1": [100]})]))
    _expect_panel_422({"plasmid_length": 1000})
