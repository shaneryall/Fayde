# fayde-detect

The `fayde-detect` crate finds stars and streaks in images for the
[Fayde](https://github.com/shaneryall/Fayde) astrometric plate solver. It measures the width of
the stars, the background, and the noise, so you don't set any of them. For each star, it
returns a position, a position error, a flux, and the probability that the star is real. It
sets aside hot pixels and cosmic-ray hits, and it traces satellite streaks and reports their
ends, angle, width, and flux.

The crate handles images of point sources, such as images taken with sidereal tracking, and
images whose stars are straight trails of one length and direction, such as images taken while
the telescope tracks a satellite. For trailed stars, it reports the trail and each star's
position error along and across it.

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
