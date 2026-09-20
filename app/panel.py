"""核心算法：覆盖全部变异型的最小酶对面板（精确最小集合覆盖）。

所有函数均为纯函数，不依赖任何外部服务或全局状态。

覆盖规则
--------
对每种酶对 (a, b) 与每个变异 v，复用 :mod:`app.digest` 的圆环双酶切与
条带聚类规则：若该酶对在"参考型 vs 变异 v"间产生至少一条可区分条带，
则称该酶对**覆盖**变异 v。任一样本有效切点不足 2 个的（酶对, 变异）
组合不可用，不产生覆盖；酶未提供某变异的切点时按空列表处理。

最小面板
--------
在所有覆盖全部变异的酶对集合中精确求出规模最小者（不使用贪心近似）：
按集合规模 k = 1, 2, … 递增做迭代加深搜索。最小面板中每种酶对在被
选中时必然带来至少一个此前未覆盖的变异（否则移除它仍全覆盖，与最小性
矛盾），据此剪枝不会丢失任何最小解。同规模方案依次按可区分条带总数
降序、排序后的 pair_key 列表字典序升序选定。
"""
from itertools import combinations

from app.digest import _pair_key, build_differences, cluster_bands, digest_sample
from app.models import (
    DifferenceOut,
    PanelPairOut,
    PanelRequest,
    PanelResponse,
    PanelVariantCoverageOut,
)


class PairCoverage:
    """一种酶对的覆盖情况：覆盖的变异位掩码、可区分条带总数与逐变异差异依据。"""

    __slots__ = ("pair", "key", "mask", "bands", "differences")

    def __init__(self, pair: list[str], key: str, mask: int, bands: int,
                 differences: dict[str, list[DifferenceOut]]):
        self.pair = pair            # 排序后的酶名对
        self.key = key              # pair_key
        self.mask = mask            # 覆盖变异的位掩码（第 i 位对应 variants[i]）
        self.bands = bands          # 覆盖的全部变异上的可区分条带总数
        self.differences = differences  # 变异名 -> 该变异上的差异条带依据


# --------------------------------------------------------------- 覆盖计算

def _pair_variant_differences(req: PanelRequest, a: str, b: str, variant: str,
                              ref_sites: dict[str, list[int]],
                              var_sites: dict[str, dict[str, list[int]]],
                              ) -> list[DifferenceOut] | None:
    """某酶对在某变异上的差异条带；切点不足或无差异（不覆盖）时返回 None。"""
    site_map = {
        enzyme: {
            "reference": ref_sites[enzyme],
            "variant": var_sites[enzyme].get(variant, []),
        }
        for enzyme in (a, b)
    }
    samples = {
        sample: digest_sample(
            site_map, a, b, sample, req.plasmid_length, req.min_visible_fragment,
        )
        for sample in ("reference", "variant")
    }
    if len(samples["reference"].cuts) < 2 or len(samples["variant"].cuts) < 2:
        return None
    bands = cluster_bands(samples, req.band_resolution_tolerance, req.min_visible_fragment)
    return build_differences(bands) or None


def compute_pair_coverage(req: PanelRequest, variants: list[str]) -> list[PairCoverage]:
    """枚举所有酶对，计算每个酶对覆盖的变异集合与逐变异差异依据。"""
    vindex = {v: i for i, v in enumerate(variants)}
    ref_sites = {e.name: sorted(e.reference_sites) for e in req.enzymes}
    var_sites = {
        e.name: {v: sorted(sites) for v, sites in e.variant_sites.items()}
        for e in req.enzymes
    }

    candidates: list[PairCoverage] = []
    for a, b in combinations([e.name for e in req.enzymes], 2):
        pair, key = _pair_key(a, b)
        mask, total_bands = 0, 0
        differences: dict[str, list[DifferenceOut]] = {}
        for variant in variants:
            diffs = _pair_variant_differences(req, a, b, variant, ref_sites, var_sites)
            if diffs:
                mask |= 1 << vindex[variant]
                total_bands += len(diffs)
                differences[variant] = diffs
        if mask:
            candidates.append(PairCoverage(pair, key, mask, total_bands, differences))
    return candidates


# --------------------------------------------------------------- 精确最小面板

def _dedupe_candidates(candidates: list[PairCoverage]) -> list[PairCoverage]:
    """覆盖集合完全相同的酶对只保留最优者：条带总数降序，并列取 pair_key 小者。"""
    best: dict[int, PairCoverage] = {}
    for cand in candidates:
        if cand.mask == 0:
            continue
        cur = best.get(cand.mask)
        if cur is None or cand.bands > cur.bands or (
            cand.bands == cur.bands and cand.key < cur.key
        ):
            best[cand.mask] = cand
    return list(best.values())


def select_minimum_panel(candidates: list[PairCoverage],
                         variant_count: int) -> list[PairCoverage] | None:
    """精确求出覆盖全部变异的最小酶对面板；无法完整覆盖时返回 None。

    依次尝试规模 k = 1, 2, …；同规模方案按可区分条带总数降序、
    排序后的 pair_key 列表字典序升序选定。
    """
    if variant_count == 0:
        return []
    full = (1 << variant_count) - 1
    cands = _dedupe_candidates(candidates)
    # 覆盖变异数降序 → 条带总数降序 → pair_key 升序，提升剪枝效率
    cands.sort(key=lambda c: (-bin(c.mask).count("1"), -c.bands, c.key))
    n = len(cands)
    pops = [bin(c.mask).count("1") for c in cands]
    suffix_or = [0] * (n + 1)
    suffix_max_bands = [0] * (n + 1)
    for i in range(n - 1, -1, -1):
        suffix_or[i] = suffix_or[i + 1] | cands[i].mask
        suffix_max_bands[i] = max(suffix_max_bands[i + 1], cands[i].bands)
    if suffix_or[0] != full:
        return None

    best_keys: list[str] | None = None
    best_bands = -1
    best_chosen: list[PairCoverage] | None = None
    chosen: list[PairCoverage] = []

    def dfs(start: int, union: int, bands_sum: int, size: int) -> None:
        nonlocal best_keys, best_bands, best_chosen
        if len(chosen) == size:
            if union == full:
                keys = sorted(c.key for c in chosen)
                if bands_sum > best_bands or (
                    bands_sum == best_bands and (best_keys is None or keys < best_keys)
                ):
                    best_keys, best_bands = keys, bands_sum
                    best_chosen = list(chosen)
            return
        remaining = size - len(chosen)
        if remaining > n - start:
            return
        if (union | suffix_or[start]) != full:
            return  # 剩余候选全部选上也无法覆盖全部变异
        if bin(union).count("1") + remaining * pops[start] < variant_count:
            return  # 覆盖变异数的上界不足
        if best_chosen is not None and bands_sum + remaining * suffix_max_bands[start] < best_bands:
            return  # 条带总数上界低于已知最优
        for i in range(start, n - remaining + 1):
            cand = cands[i]
            if not cand.mask & ~union:
                continue  # 最小面板中每种酶对必须带来新覆盖的变异
            chosen.append(cand)
            dfs(i + 1, union | cand.mask, bands_sum + cand.bands, size)
            chosen.pop()

    for size in range(1, n + 1):
        dfs(0, 0, 0, size)
        if best_chosen is not None:
            return best_chosen
    return None


# --------------------------------------------------------------- 面板组装

def build_panel(req: PanelRequest) -> PanelResponse:
    variants = sorted({v for enzyme in req.enzymes for v in enzyme.variant_sites})
    candidates = compute_pair_coverage(req, variants)

    covered = 0
    for cand in candidates:
        covered |= cand.mask
    uncovered = [v for i, v in enumerate(variants) if not (covered >> i) & 1]

    base = dict(
        plasmid_length=req.plasmid_length,
        min_visible_fragment=req.min_visible_fragment,
        band_resolution_tolerance=req.band_resolution_tolerance,
        variants=variants,
        total_enzyme_pairs=len(req.enzymes) * (len(req.enzymes) - 1) // 2,
    )
    if uncovered:
        # 任一变异无可用酶对：不返回部分面板
        return PanelResponse(
            **base, selected_pairs=[], variant_coverage=[], uncovered_variants=uncovered,
        )

    chosen = select_minimum_panel(candidates, len(variants))
    selected = sorted(chosen, key=lambda c: c.key)
    vindex = {v: i for i, v in enumerate(variants)}
    selected_pairs = [
        PanelPairOut(
            enzyme_pair=cand.pair,
            pair_key=cand.key,
            covers_variants=[v for v in variants if (cand.mask >> vindex[v]) & 1],
            distinguishable_band_count=cand.bands,
        )
        for cand in selected
    ]

    coverage: list[PanelVariantCoverageOut] = []
    for variant in variants:
        # 选中酶对中该变异可区分条带最多者作为依据（并列取 pair_key 小者；
        # selected 已按 pair_key 升序，max 保留首个最大值）
        best = max(
            (cand for cand in selected if (cand.mask >> vindex[variant]) & 1),
            key=lambda cand: len(cand.differences[variant]),
        )
        diffs = best.differences[variant]
        coverage.append(PanelVariantCoverageOut(
            variant=variant,
            enzyme_pair=best.pair,
            pair_key=best.key,
            distinguishable_band_count=len(diffs),
            differences=diffs,
        ))
    return PanelResponse(
        **base,
        selected_pairs=selected_pairs,
        variant_coverage=coverage,
        uncovered_variants=[],
    )
