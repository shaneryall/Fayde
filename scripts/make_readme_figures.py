# /// script
# requires-python = ">=3.10"
# dependencies = ["numpy", "pillow>=10.1"]
# ///
r"""Animated README figures: `fayde detect` on a synthetic frame, its report drawn over the image.

    uv run scripts/make_readme_figures.py points --fayde FAYDE \
        -o docs/images/detect-point-stars.gif
    uv run scripts/make_readme_figures.py trailed --fayde FAYDE \
        -o docs/images/detect-trailed-stars.gif

FAYDE is the `fayde` command to run, for example `target/release/fayde`. The frame is synthetic,
generated from random seed 7:

- points: 1024 × 1024 pixels, 700 Gaussian stars (σ 1.6 px), three satellite streaks, 120
  cosmic-ray hits, a sloped sky of 200 counts, Poisson noise and 6 counts of read noise;
- trailed: the same sky and noise, 350 stars trailed 48 px at 32° as when the telescope tracks a
  satellite, the satellite itself as a point source, and 120 cosmic-ray hits.

Every marker, count and measurement drawn over the image comes from the JSON report that FAYDE
prints; none comes from the simulation. The captions also describe the simulated scene, and the
subtitle gives the frame size and the version that `FAYDE --version` prints. The animation of a
streak being traced is drawing only.
"""

import argparse
import json
import os
import subprocess
import tempfile

import numpy as np
from PIL import Image, ImageDraw, ImageFont

N = 1024
SKY, READ_NOISE = 200.0, 6.0
PANEL, BAR = 760, 64


def frame_points() -> np.ndarray:
    rng = np.random.default_rng(7)
    yy, xx = np.mgrid[0:N, 0:N].astype(np.float64)
    img = SKY + 0.03 * xx - 0.02 * yy
    for _ in range(700):
        flux = 300 * 400 ** rng.uniform() ** 1.6
        spot(img, xx, yy, rng.uniform(5, N - 5), rng.uniform(5, N - 5), flux, 1.6)
    noise = np.sqrt(SKY + READ_NOISE**2)
    for x0, y0, x1, y1, snr in (
        (-30, 180, 1060, 760, 22),
        (820, -20, 600, 1050, 9),
        (210, 820, 470, 905, 14),
    ):
        line(img, xx, yy, x0, y0, x1, y1, snr * noise, 1.8, 1.5, window=False)
    img = rng.poisson(np.clip(img, 0, None)).astype(np.float64) + rng.normal(
        0, READ_NOISE, img.shape
    )
    cosmic_rays(img, rng)
    return img


def frame_trailed() -> np.ndarray:
    rng = np.random.default_rng(7)
    yy, xx = np.mgrid[0:N, 0:N].astype(np.float64)
    img = SKY + 0.03 * xx - 0.02 * yy
    length, angle = 48.0, np.deg2rad(32.0)
    hx, hy = 0.5 * length * np.cos(angle), 0.5 * length * np.sin(angle)
    for _ in range(350):
        flux = 900 * 300 ** rng.uniform() ** 1.6
        x, y = rng.uniform(-20, N + 20), rng.uniform(-20, N + 20)
        line(img, xx, yy, x - hx, y - hy, x + hx, y + hy, flux / length, 1.6, 1.2, window=True)
    spot(img, xx, yy, 540.0, 470.0, 60000.0, 1.6)
    img = rng.poisson(np.clip(img, 0, None)).astype(np.float64) + rng.normal(
        0, READ_NOISE, img.shape
    )
    cosmic_rays(img, rng)
    return img


def spot(img, xx, yy, x, y, flux, s):
    r = int(4 * s + 3)
    x0, x1 = max(int(x) - r, 0), min(int(x) + r + 1, N)
    y0, y1 = max(int(y) - r, 0), min(int(y) + r + 1, N)
    gx, gy = xx[y0:y1, x0:x1] - x, yy[y0:y1, x0:x1] - y
    img[y0:y1, x0:x1] += flux / (2 * np.pi * s * s) * np.exp(-(gx * gx + gy * gy) / (2 * s * s))


def line(img, xx, yy, x0, y0, x1, y1, per_px, s, soft, window):
    """A segment with soft ends convolved with a Gaussian of σ `s`, `per_px` counts per pixel."""
    length = np.hypot(x1 - x0, y1 - y0)
    ux, uy = (x1 - x0) / length, (y1 - y0) / length
    if window:
        r = int(length / 2 + 6 * s)
        cx, cy = (x0 + x1) / 2, (y0 + y1) / 2
        xa, xb = max(int(cx) - r, 0), min(int(cx) + r + 1, N)
        ya, yb = max(int(cy) - r, 0), min(int(cy) + r + 1, N)
        if xa >= xb or ya >= yb:
            return
    else:
        xa, xb, ya, yb = 0, N, 0, N
    X, Y = xx[ya:yb, xa:xb], yy[ya:yb, xa:xb]
    t = (X - x0) * ux + (Y - y0) * uy
    d = -(X - x0) * uy + (Y - y0) * ux
    along = 0.5 * (np.tanh(t / soft) - np.tanh((t - length) / soft))
    img[ya:yb, xa:xb] += per_px / (np.sqrt(2 * np.pi) * s) * np.exp(-d * d / (2 * s * s)) * along


def cosmic_rays(img, rng):
    for _ in range(120):
        x, y = rng.integers(2, N - 2, 2)
        img[y, x] += rng.uniform(800, 6000)
        if rng.uniform() < 0.4:
            ang = rng.uniform(0, np.pi)
            for k in range(1, rng.integers(2, 6)):
                xi, yi = int(x + k * np.cos(ang)), int(y + k * np.sin(ang))
                if 0 <= xi < N and 0 <= yi < N:
                    img[yi, xi] += rng.uniform(600, 4000)


def write_fits(path, img):
    cards = [
        "SIMPLE  =                    T",
        "BITPIX  =                  -32",
        "NAXIS   =                    2",
        f"NAXIS1  = {N:20d}",
        f"NAXIS2  = {N:20d}",
        "END",
    ]
    header = "".join(c.ljust(80) for c in cards).ljust(2880).encode()
    data = img.astype(">f4").tobytes()
    with open(path, "wb") as f:
        f.write(header + data + b"\0" * (-len(data) % 2880))


class Canvas:
    """Frames of the animation: the stretched image, the overlay, and a caption bar."""

    def __init__(self, img, subtitle):
        lo, hi = np.percentile(img, [5, 99.7])
        v = np.arcsinh(10 * np.clip((img - lo) / (hi - lo), 0, 1)) / np.arcsinh(10)
        # FITS row 1 is the bottom row: flip so the image reads as in a FITS viewer.
        self.base = (
            Image.fromarray((v[::-1] * 255).astype(np.uint8))
            .resize((PANEL, PANEL), Image.LANCZOS)
            .convert("RGB")
        )
        self.font = ImageFont.load_default(size=20)
        self.small = ImageFont.load_default(size=15)
        self.subtitle = subtitle
        self.frames, self.durations = [], []

    def add(self, draw, caption, ms):
        panel = self.base.copy()
        draw(ImageDraw.Draw(panel, "RGBA"))
        im = Image.new("RGB", (PANEL, PANEL + BAR), (12, 12, 16))
        im.paste(panel, (0, 0))
        d = ImageDraw.Draw(im)
        font = self.font
        for size in (18, 17, 16):
            if d.textlength(caption, font=font) <= PANEL - 28:
                break
            font = ImageFont.load_default(size=size)
        d.text((14, PANEL + 10), caption, font=font, fill=(235, 235, 235))
        d.text((14, PANEL + 38), self.subtitle, font=self.small, fill=(140, 140, 150))
        self.frames.append(im)
        self.durations.append(ms)

    def save(self, path):
        # One 128-color palette, taken from the first, middle and last frames together, keeps the
        # file under 1 MB without shifting the colors of the image or the overlay.
        picks = [self.frames[0], self.frames[len(self.frames) // 2], self.frames[-1]]
        sample = Image.new("RGB", (PANEL * len(picks), PANEL + BAR))
        for i, f in enumerate(picks):
            sample.paste(f, (PANEL * i, 0))
        palette = sample.quantize(colors=128, method=Image.Quantize.FASTOCTREE)
        frames = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in self.frames]
        frames[0].save(
            path,
            save_all=True,
            append_images=frames[1:],
            duration=self.durations,
            loop=0,
            optimize=True,
        )


def pos(x, y):
    """FITS 1-based pixel to panel coordinates."""
    s = PANEL / N
    return np.array(((x - 0.5) * s, (N - (y - 0.5)) * s))


def cross(d, p, color, k=4):
    d.line([tuple(p - k), tuple(p + k)], fill=color, width=2)
    d.line([(p[0] - k, p[1] + k), (p[0] + k, p[1] - k)], fill=color, width=2)


def dot(d, p, color, r=1.5):
    d.ellipse([p[0] - r, p[1] - r, p[0] + r, p[1] + r], fill=color)


def ring(d, p, r, color, width=1):
    d.ellipse([p[0] - r, p[1] - r, p[0] + r, p[1] + r], outline=color, width=width)


def label(d, canvas, at, text, color):
    tw = d.textlength(text, font=canvas.small)
    x0 = min(max(at[0] - tw / 2, 4), PANEL - tw - 12)
    y0 = min(max(at[1] - 10, 4), PANEL - 26)
    d.rounded_rectangle(
        [x0 - 6, y0 - 3, x0 + tw + 6, y0 + 20], 5, fill=(0, 0, 0, 170), outline=color
    )
    d.text((x0, y0), text, font=canvas.small, fill=color)


STREAK_COLORS = [(255, 196, 0), (0, 220, 255), (255, 90, 200), (160, 255, 90)]
SHARP, OTHER, STAR = (255, 60, 60), (255, 140, 0), (80, 255, 140)
CANDIDATE = (160, 255, 90)
ORANGE_LABELS = {"streak": "on a streak", "extended": "extended"}
WORDS = ("zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine")


def num(n):
    """A count as running text writes it: words below 10, numerals from 10 up."""
    return WORDS[n] if n < 10 else str(n)


def plural(n, one, many=None):
    return f"{num(n)} {one if n == 1 else (many or one + 's')}"


def set_aside(sharp, on_streak):
    # A count below 10 stays a numeral in a sentence that has one of 10 or more.
    say = str if max(len(sharp), len(on_streak)) >= 10 else num
    noun = "detection" if len(sharp) == 1 else "detections"
    text = f"Fayde sets aside {say(len(sharp))} {noun} as too sharp (red)"
    if on_streak:
        text += f" and {say(len(on_streak))} on a streak (orange)"
    return text + "."


def orange(*groups):
    """What the orange marks are, from the classes in the report."""
    classes = sorted({s["class"] for g in groups for s in g}, key=lambda c: (c != "streak", c))
    return " or ".join(ORANGE_LABELS.get(c, c) for c in classes)


def key(*entries):
    """The color key for the last frame, without colors that mark nothing."""
    return " ".join(f"{color}: {what}." for color, what in entries if what)


def draw_streak(d, canvas, i, st, frac=1.0, show_label=True, color=None):
    c = color or STREAK_COLORS[i % len(STREAK_COLORS)]
    a, b = pos(st["x0"], st["y0"]), pos(st["x1"], st["y1"])
    e = a + frac * (b - a)
    d.line([tuple(a), tuple(e)], fill=c + (90,), width=9)
    d.line([tuple(a), tuple(e)], fill=c + (255,), width=2)
    for p in (a,) + ((b,) if frac >= 1 else ()):
        ring(d, p, 5, c, 2)
    if show_label and frac >= 1:
        m = a + (0.22, 0.78, 0.5, 0.35)[i % 4] * (b - a)
        n = np.array([-(b - a)[1], (b - a)[0]]) / np.hypot(*(b - a))
        label(
            d,
            canvas,
            m + 22 * (n if n[1] < 0 else -n),
            f"#{i + 1}  {st['length_px']:.0f} px  SNR {st['snr']:.0f}",
            c,
        )


def figure_points(report, img, out, subtitle):
    canvas = Canvas(img, subtitle)
    stars = [s for s in report["sources"] if s["class"] == "star"]
    others = [s for s in report["sources"] if s["class"] != "star"]
    sharp = [r for r in report["rejected"] if r["class"] == "sharp"]
    on_streak = [r for r in report["rejected"] if r["class"] == "streak"]
    rest = [r for r in report["rejected"] if r["class"] not in ("sharp", "streak")]
    streaks = report["streaks"]

    def marks(d, a=255):
        for r in sharp:
            cross(d, pos(r["x"], r["y"]), SHARP + (a,))
        for r in on_streak:
            dot(d, pos(r["x"], r["y"]), OTHER + (a,))
        for r in rest + others:
            ring(d, pos(r["x"], r["y"]), 7, OTHER + (a,), 2)

    def star_rings(d, a=150):
        for s in stars:
            ring(d, pos(s["x"], s["y"]), 5, STAR + (a,))

    canvas.add(lambda d: None, "The frame before detection.", 1800)
    canvas.add(marks, set_aside(sharp, on_streak), 2000)
    kept = "a star" if len(stars) == 1 else "stars"
    canvas.add(
        lambda d: (marks(d, 110), star_rings(d)),
        f"Fayde keeps {plural(len(stars), 'detection')} as {kept} (green).",
        1700,
    )
    done = []
    for i, st in enumerate(streaks):
        for k in range(1, 15):
            canvas.add(
                lambda d, i=i, k=k, st=st: (
                    marks(d, 60),
                    star_rings(d, 60),
                    [draw_streak(d, canvas, j, streaks[j]) for j in done],
                    draw_streak(d, canvas, i, st, k / 14, False),
                ),
                f"Streak {i + 1} of {len(streaks)}",
                45,
            )
        done.append(i)
        canvas.add(
            lambda d, done=tuple(done): (
                marks(d, 60),
                star_rings(d, 60),
                [draw_streak(d, canvas, j, streaks[j]) for j in done],
            ),
            f"Streak {i + 1} is {st['length_px']:.0f} px long at {st['angle_deg']:.0f}°, "
            f"{st['width_sigma_px']:.1f} px wide (Gaussian sigma).",
            1300,
        )
    canvas.add(
        lambda d: (
            marks(d, 110),
            star_rings(d, 110),
            [draw_streak(d, canvas, j, streaks[j]) for j in done],
        ),
        key(
            ("Green", "stars" if stars else ""),
            ("Red", "too sharp" if sharp else ""),
            ("Orange", orange(on_streak, rest, others)),
        ),
        3500,
    )
    canvas.save(out)


def figure_trailed(report, img, out, subtitle):
    canvas = Canvas(img, subtitle)
    trail = report["trail"]
    length, ang = trail["length_px"], np.deg2rad(trail["angle_deg"])
    u = np.array((np.cos(ang), -np.sin(ang))) * PANEL / N
    stars = sorted((s for s in report["sources"] if s["class"] == "star"), key=lambda s: -s["flux"])
    others = [s for s in report["sources"] if s["class"] != "star"]
    sharp = [r for r in report["rejected"] if r["class"] == "sharp"]
    points = [r for r in report["rejected"] if r["class"] == "point"]
    on_streak = [r for r in report["rejected"] if r["class"] == "streak"]
    rest = [r for r in report["rejected"] if r["class"] not in ("sharp", "point", "streak")]
    streaks = report["streaks"]

    def marks(d, a=255):
        for r in sharp:
            cross(d, pos(r["x"], r["y"]), SHARP + (a,))
        for r in on_streak + rest + others:
            ring(d, pos(r["x"], r["y"]), 7, OTHER + (a,), 2)

    def trails(d, n, a=255):
        for s in stars[:n]:
            c = pos(s["x"], s["y"])
            p, q = c - 0.5 * length * u, c + 0.5 * length * u
            d.line([tuple(p), tuple(q)], fill=(0, 220, 255, a // 6), width=6)
            d.line([tuple(p), tuple(q)], fill=(0, 220, 255, a), width=1)

    def point_sources(d):
        for r in points:
            p = pos(r["x"], r["y"])
            ring(d, p, 12, (255, 196, 0), 3)
            label(d, canvas, p + np.array((0, -26)), f"point  SNR {r['snr']:.0f}", (255, 196, 0))

    def candidates(d):
        for i, st in enumerate(streaks):
            draw_streak(d, canvas, i, st, color=CANDIDATE)

    canvas.add(
        lambda d: None, "The telescope tracks a satellite, so every star is a short trail.", 2000
    )
    canvas.add(marks, set_aside(sharp, on_streak), 2000)
    for k in range(1, 13):
        n = round(len(stars) * (k / 12) ** 1.5)
        canvas.add(
            lambda d, n=n: (marks(d, 70), trails(d, n)),
            f"Trails, brightest first: {n} of {len(stars)}",
            70,
        )
    canvas.add(
        lambda d: (marks(d, 70), trails(d, len(stars))),
        f"Fayde measures a position for each of the {num(len(stars))} trails it finds.",
        1600,
    )
    canvas.add(
        lambda d: (marks(d, 70), trails(d, len(stars))),
        f"A fit to {num(trail['fitted'])} trails gives the whole frame one length and angle: "
        f"{length:.1f} px at {trail['angle_deg']:.1f}°.",
        2400,
    )
    canvas.add(
        lambda d: (trails(d, len(stars), 90), point_sources(d)),
        f"Fayde reports {plural(len(points), 'detection')} that "
        f"{'doesn' if len(points) == 1 else 'don'}'t trail as "
        f"{'a point source' if len(points) == 1 else 'point sources'} (yellow).",
        2200,
    )
    if streaks:
        canvas.add(
            lambda d: (trails(d, len(stars), 90), point_sources(d), candidates(d)),
            f"Fayde also reports {plural(len(streaks), 'streak candidate')} (green).",
            2200,
        )
    canvas.add(
        lambda d: (marks(d, 110), trails(d, len(stars), 160), point_sources(d), candidates(d)),
        key(
            ("Cyan", "trails"),
            ("Yellow", "point sources" if points else ""),
            (
                "Green",
                "" if not streaks else "streak candidates" if streaks[1:] else "streak candidate",
            ),
            ("Red", "too sharp" if sharp else ""),
            ("Orange", orange(on_streak, rest, others)),
        ),
        3500,
    )
    canvas.save(out)


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    ap.add_argument("kind", choices=["points", "trailed"])
    ap.add_argument("--fayde", required=True, help="the fayde command to run")
    ap.add_argument("-o", "--output", required=True, help="GIF to write")
    ap.add_argument("--report", help="also keep the JSON report here")
    a = ap.parse_args()
    img = frame_points() if a.kind == "points" else frame_trailed()
    with tempfile.TemporaryDirectory() as tmp:
        fits = os.path.join(tmp, "frame.fits")
        write_fits(fits, img)
        text = subprocess.run(
            [a.fayde, "detect", "--image", fits], check=True, capture_output=True, text=True
        ).stdout
    report = json.loads(text)
    version = subprocess.run(
        [a.fayde, "--version"], check=True, capture_output=True, text=True
    ).stdout.split()[1]
    subtitle = (
        f"Synthetic frame of 1024 x 1024 pixels, run through fayde detect (version {version})"
    )
    if a.report:
        with open(a.report, "w") as f:
            f.write(text)
    (figure_points if a.kind == "points" else figure_trailed)(report, img, a.output, subtitle)


if __name__ == "__main__":
    main()
