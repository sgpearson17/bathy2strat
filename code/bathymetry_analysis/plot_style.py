"""Shared plot styling utilities."""

from __future__ import annotations

import matplotlib as mpl
from matplotlib.font_manager import FontProperties

_FONT_FAMILY = "Arial"
_FONT_WEIGHT = "bold"
_FONT_STYLE = "italic"


def get_plot_font() -> FontProperties:
    return FontProperties(family=_FONT_FAMILY, weight=_FONT_WEIGHT, style=_FONT_STYLE)


def apply_axes_font(ax, font: FontProperties | None = None) -> None:
    font = font or get_plot_font()
    ax.title.set_fontproperties(font)
    ax.xaxis.label.set_fontproperties(font)
    ax.yaxis.label.set_fontproperties(font)
    for label in ax.get_xticklabels() + ax.get_yticklabels():
        label.set_fontproperties(font)
    for text in ax.texts:
        text.set_fontproperties(font)


def apply_global_style() -> None:
    mpl.rcParams["font.family"] = _FONT_FAMILY
    mpl.rcParams["font.weight"] = _FONT_WEIGHT
    mpl.rcParams["font.style"] = _FONT_STYLE
    mpl.rcParams["axes.titleweight"] = _FONT_WEIGHT
    mpl.rcParams["axes.labelweight"] = _FONT_WEIGHT
