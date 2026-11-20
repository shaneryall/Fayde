# fayde-api

The `fayde-api` crate builds requests and reports for the
[Fayde](https://github.com/shaneryall/Fayde) astrometric plate solver. The command-line tool and
Python package share this crate. It builds a solve request from your options, an instrument
profile, and a frame header, in that order of precedence. It then runs the request and returns
the JSON report.

For more information about Fayde and how to install the `fayde` command or Python package, see
the [Fayde README](https://github.com/shaneryall/Fayde#readme).

All Fayde crates share one release version. Each crate requires exactly that version of the
other Fayde crates.

For more information about the algorithms, see the
[references](https://github.com/shaneryall/Fayde/blob/main/docs/references.md).

## License

Fayde uses the GNU General Public License, version 3 or later. For license terms and
attribution, see the `LICENSE` and `NOTICE` files in this crate. If you use Fayde in a publication, cite it as described in
[`CITATION.cff`](https://github.com/shaneryall/Fayde/blob/main/CITATION.cff).
