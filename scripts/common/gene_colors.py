#!/usr/bin/env python3
"""
gene_colors.py
Deterministic randomised colour assignment per gene name, replacing a
hand-maintained hex-code dict (which doesn't scale -- every new gene name
added to the pipeline needed a manually-picked colour).

Colours are generated in HSL space with a minimum-distance rejection
threshold: each new colour is retried until it's at least MIN_DISTANCE away
(Euclidean, in RGB 0-255 space) from every colour already assigned -- i.e.
it can't fall "within a circle" of an existing colour. Deterministic (fixed
seed + genes sorted before assignment) so the same gene always gets the
same colour across separate runs/scripts, as long as the same gene set is
passed in.
"""

import colorsys
import random

SEED = 42
MIN_DISTANCE = 90       # Euclidean RGB distance an new colour must clear
MAX_ATTEMPTS = 500       # per colour, before relaxing the threshold
SATURATION_RANGE = (0.55, 0.80)
LIGHTNESS_RANGE = (0.40, 0.60)


def _rgb_distance(a: tuple[int, int, int], b: tuple[int, int, int]) -> float:
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def _hsl_to_hex(h: float, s: float, l: float) -> tuple[str, tuple[int, int, int]]:
    r, g, b = colorsys.hls_to_rgb(h / 360, l, s)
    rgb = (round(r * 255), round(g * 255), round(b * 255))
    return "#{:02X}{:02X}{:02X}".format(*rgb), rgb


DARKNESS_THRESHOLD = 128  # perceived-luminance cutoff (0-255); below this, label text goes white


def is_dark(hex_code: str, threshold: float = DARKNESS_THRESHOLD) -> bool:
    """True if hex_code is dark enough that black label text would be hard
    to read on it (perceived luminance, ITU-R BT.601 weights)."""
    hex_code = hex_code.lstrip("#")
    r, g, b = (int(hex_code[i:i + 2], 16) for i in (0, 2, 4))
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return luminance < threshold


def text_color_for(hex_code: str) -> str:
    """Label colour to use on top of hex_code's fill -- white if hex_code is
    dark, black otherwise."""
    return "white" if is_dark(hex_code) else "black"


def generate_gene_colors(gene_names, seed: int = SEED, min_distance: float = MIN_DISTANCE) -> dict[str, str]:
    """gene_names: any iterable of gene name strings (case-insensitive).
    Returns {gene_name.lower(): "#RRGGBB"}, sorted-then-assigned so the
    mapping is stable across runs regardless of input order."""
    rng = random.Random(seed)
    used_rgb: list[tuple[int, int, int]] = []
    colors: dict[str, str] = {}

    for gene in sorted({g.lower() for g in gene_names}):
        threshold = min_distance
        chosen_hex = chosen_rgb = None
        # Relax the threshold in steps if we can't find a far-enough colour
        # in time (happens once the palette gets crowded) -- better to
        # place a slightly-too-close colour than loop forever.
        while chosen_hex is None:
            for _ in range(MAX_ATTEMPTS):
                h = rng.uniform(0, 360)
                s = rng.uniform(*SATURATION_RANGE)
                l = rng.uniform(*LIGHTNESS_RANGE)
                hex_code, rgb = _hsl_to_hex(h, s, l)
                if all(_rgb_distance(rgb, u) >= threshold for u in used_rgb):
                    chosen_hex, chosen_rgb = hex_code, rgb
                    break
            if chosen_hex is None:
                threshold *= 0.85  # relax and retry
                if threshold < 10:
                    # give up relaxing further, just take the last attempt
                    chosen_hex, chosen_rgb = hex_code, rgb

        used_rgb.append(chosen_rgb)
        colors[gene] = chosen_hex

    return colors


if __name__ == "__main__":
    demo_genes = ["flga", "flgb", "flgc", "flik", "mota", "motb", "pfla", "maf_2356", "maf7"]
    for gene, hexcode in generate_gene_colors(demo_genes).items():
        print(f"{gene:10s} {hexcode}")
