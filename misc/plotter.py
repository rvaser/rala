#!/usr/bin/env python

from __future__ import print_function
import os, sys, argparse, json, matplotlib.pyplot

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#*******************************************************************************

class Plotter:
    def __init__(self, knots, out_directory):
        self.knots = knots
        self.out_directory = out_directory + "/"
        if (not os.path.isdir(self.out_directory)):
            eprint("[rala::Plotter::__init__] error: invalid out directory!")
            sys.exit(1)

    def __enter__(self):
        pass

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    @staticmethod
    def plot_pile(pile, orientation, overlap_length, type, title, ax):
        if ("y" not in pile or "~y" not in pile or "b" not in pile or
            "e" not in pile or "h" not in pile or "m" not in pile or
            "p10" not in pile):
            eprint("[rala::Plotter::plot_pile] error: incomplete pile!")
            sys.exit(1);

        x = xrange(len(pile["y"]))

        if (orientation == 0):
            ax.plot(x, pile["y"], label="y")
            if (pile["~y"]):
                ax.plot(x, pile["~y"], label="~y")
            for slope in pile["h"]:
                ax.axvline(slope, color="r", linestyle=":")
            begin = pile["b"]
            end = pile["e"]

        else:
            ax.plot(x, list(reversed(pile["y"])), label="y")
            if (pile["~y"]):
                ax.plot(x, list(reversed(pile["~y"])), label="~y")
            for slope in pile["h"]:
                ax.axvline(len(x) - slope, color="r", linestyle=":")
            begin = len(x) - pile["e"]
            end = len(x) - pile["b"]

        if (type == "p"):
            ax.axvline(begin + overlap_length, label="o", color="k", linestyle="--")
        else:
            ax.axvline(end - overlap_length, label="o", color="k", linestyle="--")

        ax.axhline(int(pile["m"]), label="m", color="m", linestyle=":")
        ax.axhline(int(pile["p10"]), label="p10", color="c", linestyle=":")
        ax.set_title(title)

    def run(self):
        try:
            k = open(self.knots)
        except Exception:
            eprint("[rala::Plotter::run] error: unable to open file {}!".format(self.knots))
            sys.exit(1)

        try:
            data = json.load(k)
        except Exception:
            eprint("[rala::Plotter::run] error: file is not in JSON format!")
            sys.exit(1)

        if ("knots" not in data or not data["knots"] or\
            "piles" not in data or not data["piles"]):
            eprint("[rala::Plotter::run] error: incomplete input file!")
            sys.exit(1)

        for knot in data["knots"]:
            if (knot not in data["piles"]):
                eprint("[rala::Plotter::run] error: missing pile {}!".format(knot))
                sys.exit(1)

            if ("n" not in data["knots"][knot] or "p" not in data["knots"][knot] or\
                "s" not in data["knots"][knot]):
                eprint("[rala::Plotter::run] error: incomplete knot!")
                sys.exit(1)

            num_plots = len(data["knots"][knot]["p"]) + len(data["knots"][knot]["s"])
            if (num_plots == 0):
                continue

            figure, axes = matplotlib.pyplot.subplots(num_plots, 2, sharey="all",\
                figsize=(15,15))

            ax_row = 0
            for prefix in data["knots"][knot]["p"]:
                if (prefix[0] not in data["piles"]):
                    eprint("[rala::Plotter::run] error: missing pile {}!".format(prefix[0]))
                    sys.exit(1)

                Plotter.plot_pile(data["piles"][prefix[0]], prefix[2], prefix[3],\
                    "s", prefix[0] + " - " + prefix[1], axes[ax_row, 0])
                Plotter.plot_pile(data["piles"][knot], 0, prefix[3],\
                    "p", knot, axes[ax_row, 1])
                ax_row += 1

            for suffix in data["knots"][knot]["s"]:
                if (suffix[0] not in data["piles"]):
                    eprint("[rala::Plotter::run] error: missing pile {}!".format(suffix[0]))
                    sys.exit(1)

                Plotter.plot_pile(data["piles"][knot], 0, suffix[3],\
                    "s", knot, axes[ax_row, 0])
                Plotter.plot_pile(data["piles"][suffix[0]], suffix[3], suffix[3],\
                    "p", suffix[0] + " - " + suffix[1], axes[ax_row, 1])
                ax_row += 1

            figure.text(0.5, 0.04, "base", ha="center")
            figure.text(0.04, 0.5, "coverage", va="center", rotation="vertical")
            matplotlib.pyplot.legend(loc="best")
            matplotlib.pyplot.savefig(self.out_directory +\
                str(data["knots"][knot]['n']) + "_knots.png")
            matplotlib.pyplot.close(figure)

        sys.exit(1);


#*******************************************************************************

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Plotter is a handy tool for
        drawing leftover knots in an assembly graph constructed with rala""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("knots", help="""input file in JSON format containing
        information about read piles and overlaps between them""")
    parser.add_argument("-o", "--out-directory", default=os.getcwd(),
        help="""path in which plotted images will be saved""")

    args = parser.parse_args()

    plotter = Plotter(args.knots, args.out_directory)

    with plotter:
        plotter.run()
