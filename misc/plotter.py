#!/usr/bin/env python

from __future__ import print_function
import os, sys, argparse, json, matplotlib.pyplot

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

#*******************************************************************************

class Plotter:
    def __init__(self, data_path, ylimit, out_path):
        self.data_path = data_path
        self.ylimit = int(ylimit)
        self.out_path = out_path + "/"
        if (not os.path.isdir(self.out_path)):
            eprint("[rala::Plotter::__init__] error: invalid out directory!")
            sys.exit(1)

    def __enter__(self):
        pass

    def __exit__(self, exception_type, exception_value, traceback):
        pass

    @staticmethod
    def plot_pile(pile, orientation, overlap_length, type, title, ylimit, ax):
        if ("y" not in pile or "b" not in pile or "e" not in pile or
            "h" not in pile or "m" not in pile or "p10" not in pile):
            eprint("[rala::Plotter::plot_pile] error: incomplete pile!")
            sys.exit(1);

        x = xrange(len(pile["y"]))
        ax.set_ylim([0, ylimit])

        if (orientation == 0):
            ax.plot(x, pile["y"], label="y")
            for slope in pile["h"]:
                ax.axvline(slope, color="r", linestyle=":")
            begin = pile["b"]
            end = pile["e"]
        else:
            ax.plot(x, list(reversed(pile["y"])), label="y")
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
            d = open(self.data_path)
        except Exception:
            eprint("[rala::Plotter::run] error: unable to open file {}!".format(self.data_path))
            sys.exit(1)

        try:
            data = json.load(d)
        except Exception:
            eprint("[rala::Plotter::run] error: file is not in JSON format!")
            sys.exit(1)

        if ("piles" not in data or not data["piles"]):
            eprint("[rala::Plotter::run] error: incomplete input file!")
            sys.exit(1)

        if ("nodes" not in data or not data["nodes"]):
            for pile in data["piles"]:
                figure, ax = matplotlib.pyplot.subplots(1, 1, figsize=(7.5, 5))

                Plotter.plot_pile(data["piles"][pile], 0, 0, "p", str(pile),\
                    self.ylimit, ax)

                figure.text(0.5, 0.04, "base", ha="center")
                figure.text(0.04, 0.5, "coverage", va="center", rotation="vertical")
                matplotlib.pyplot.legend(loc="best")
                matplotlib.pyplot.savefig(self.out_path + str(pile) + ".png")
                matplotlib.pyplot.close(figure)

            return

        for node in data["nodes"]:
            if (node not in data["piles"]):
                eprint("[rala::Plotter::run] error: missing pile {}!".format(node))
                sys.exit(1)

            if ("n" not in data["nodes"][node] or "p" not in data["nodes"][node] or\
                "s" not in data["nodes"][node]):
                eprint("[rala::Plotter::run] error: incomplete node!")
                sys.exit(1)

            num_plots = len(data["nodes"][node]["p"]) + len(data["nodes"][node]["s"])
            if (num_plots == 0):
                continue

            figure, axes = matplotlib.pyplot.subplots(num_plots, 2, sharey="all",\
                squeeze=False, figsize=(15, 5 * num_plots))

            ax_row = 0
            for prefix in data["nodes"][node]["p"]:
                if (prefix[0] not in data["piles"]):
                    eprint("[rala::Plotter::run] error: missing pile {}!".format(prefix[0]))
                    sys.exit(1)

                Plotter.plot_pile(data["piles"][prefix[0]], prefix[2], prefix[3],\
                    "s", prefix[0] + " - " + prefix[1], self.ylimit, axes[ax_row, 0])
                Plotter.plot_pile(data["piles"][node], 0, prefix[3],\
                    "p", node, self.ylimit, axes[ax_row, 1])
                ax_row += 1

            for suffix in data["nodes"][node]["s"]:
                if (suffix[0] not in data["piles"]):
                    eprint("[rala::Plotter::run] error: missing pile {}!".format(suffix[0]))
                    sys.exit(1)

                Plotter.plot_pile(data["piles"][node], 0, suffix[3],\
                    "s", node, self.ylimit, axes[ax_row, 0])
                Plotter.plot_pile(data["piles"][suffix[0]], suffix[2], suffix[3],\
                    "p", suffix[0] + " - " + suffix[1], self.ylimit, axes[ax_row, 1])
                ax_row += 1

            figure.text(0.5, 0.04, "base", ha="center")
            figure.text(0.04, 0.5, "coverage", va="center", rotation="vertical")
            matplotlib.pyplot.legend(loc="best")
            matplotlib.pyplot.savefig(self.out_path +\
                str(data["nodes"][node]['n']) + ".png")
            matplotlib.pyplot.close(figure)

#*******************************************************************************

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="""Plotter is a handy tool for
        drawing pile-ograms of nodes in an assembly graph constructed with rala""",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("data_path", help="""input file in JSON format containing
        information about read piles and overlaps between them""")
    parser.add_argument("-o", "--out-path", default=os.getcwd(),
        help="""path in which plotted images will be saved""")

    required_arguments = parser.add_argument_group('required arguments')
    required_arguments.add_argument("-y", "--ylimit", help="y axis limit",
        required=True)

    args = parser.parse_args()

    plotter = Plotter(args.data_path, args.ylimit, args.out_path)

    with plotter:
        plotter.run()
