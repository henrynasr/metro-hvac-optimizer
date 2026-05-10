"""Shared helpers: data loading, plot styling, thermal ODE slope."""

import pandas as pd
import numpy as np


def load_data(csv_path):
    """Load Paris weather CSV, parse `time` as datetime, set as index."""
    df = pd.read_csv(csv_path)
    df["time"] = pd.to_datetime(df["time"])
    df.set_index("time", inplace=True)
    return df


def style_axes(ax, title="", xlabel="", ylabel=""):
    """Apply project-standard styling to a matplotlib axis."""
    ax.set_title(title, fontsize=16)
    ax.set_xlabel(xlabel, fontsize=14)
    ax.set_ylabel(ylabel, fontsize=14)
    ax.tick_params(labelsize=11)
    ax.grid(True, alpha=0.3)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)