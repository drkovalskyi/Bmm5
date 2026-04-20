#!/usr/bin/env python3
"""Make comparison plots between current and reference NanoAOD output.

Uses PyROOT TTree::Draw to produce overlay histograms.
Reference: black, Current: magenta, line width 2.
"""

import argparse
import json
import os
import re
import sys

import ROOT
ROOT.gROOT.SetBatch(True)


def safe_filename(title):
    """Convert plot title to a safe filename."""
    name = title.lower()
    name = re.sub(r"[^a-z0-9]+", "_", name)
    return name.strip("_")


def make_plot(title, draw_expr, selection, cur_file, ref_file, output_dir):
    """Produce a single overlay plot."""
    ROOT.gStyle.SetOptStat("emr")

    c = ROOT.TCanvas("c", title, 800, 600)

    # Parse the draw expression to get unique histogram names
    # e.g., "bkmm_jpsimc_mass>>h(100,4.9,5.9)"
    draw_ref = re.sub(r">>h\(", ">>h_ref(", draw_expr)
    draw_cur = re.sub(r">>h\(", ">>h_cur(", draw_expr)

    # Reference histogram
    f_ref = ROOT.TFile.Open(ref_file)
    t_ref = f_ref.Get("Events")
    t_ref.Draw(draw_ref, selection)
    h_ref = ROOT.gDirectory.Get("h_ref").Clone("h_reference")
    h_ref.SetDirectory(0)
    f_ref.Close()

    h_ref.SetLineColor(ROOT.kBlack)
    h_ref.SetLineWidth(2)
    h_ref.SetTitle(title)

    # Current histogram
    f_cur = ROOT.TFile.Open(cur_file)
    t_cur = f_cur.Get("Events")
    t_cur.Draw(draw_cur, selection)
    h_cur = ROOT.gDirectory.Get("h_cur").Clone("h_current")
    h_cur.SetDirectory(0)
    f_cur.Close()

    h_cur.SetLineColor(ROOT.kMagenta)
    h_cur.SetLineWidth(2)

    # Y-axis range with headroom for stats boxes
    ymax = max(h_ref.GetMaximum(), h_cur.GetMaximum())
    h_ref.SetMaximum(ymax * 1.5)
    h_ref.SetMinimum(0)

    # Draw overlay
    h_ref.Draw("hist")
    h_cur.Draw("hist sames")

    # Position stats boxes
    c.Update()
    st_ref = h_ref.FindObject("stats")
    if st_ref:
        st_ref.SetTextColor(ROOT.kBlack)
        st_ref.SetLineColor(ROOT.kBlack)

    st_cur = h_cur.FindObject("stats")
    if st_cur and st_ref:
        st_cur.SetTextColor(ROOT.kMagenta)
        st_cur.SetLineColor(ROOT.kMagenta)
        # Move current stats box below reference
        height = st_ref.GetY2NDC() - st_ref.GetY1NDC()
        st_cur.SetY2NDC(st_ref.GetY1NDC())
        st_cur.SetY1NDC(st_ref.GetY1NDC() - height)
        st_cur.SetX1NDC(st_ref.GetX1NDC())
        st_cur.SetX2NDC(st_ref.GetX2NDC())

    c.Update()

    # Save
    fname = safe_filename(title) + ".png"
    fpath = os.path.join(output_dir, fname)
    c.SaveAs(fpath)
    return fname


def main():
    parser = argparse.ArgumentParser(description="Make comparison plots")
    parser.add_argument("current", help="current NanoAOD ROOT file")
    parser.add_argument("reference", help="reference NanoAOD ROOT file")
    parser.add_argument("--plots", required=True,
                        help='JSON string: {"title": ["draw_expr", "selection"], ...}')
    parser.add_argument("-o", "--output", default=".", help="output directory")
    args = parser.parse_args()

    plots = json.loads(args.plots)
    if not plots:
        return

    os.makedirs(args.output, exist_ok=True)

    for title, (draw_expr, selection) in plots.items():
        fname = make_plot(title, draw_expr, selection,
                          args.current, args.reference, args.output)
        print(f"{fname}")


if __name__ == "__main__":
    main()
