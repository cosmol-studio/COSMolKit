"""Export the editable SVG social card to the PNG used by link preview crawlers."""

import argparse
import base64
from pathlib import Path

from playwright.sync_api import sync_playwright


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--channel", help="Use an installed browser, for example msedge")
    args = parser.parse_args()
    docs_web = Path(__file__).resolve().parents[1]
    output = docs_web / "target/social-card.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    svg = (docs_web / "assets/social-card.svg").read_text(encoding="utf-8")
    logo = base64.b64encode((docs_web / "assets/logo.svg").read_bytes()).decode("ascii")
    svg = svg.replace('href="logo.svg"', f'href="data:image/svg+xml;base64,{logo}"')
    with sync_playwright() as playwright:
        browser = playwright.chromium.launch(channel=args.channel)
        page = browser.new_page(viewport={"width": 1200, "height": 630}, device_scale_factor=1)
        page.set_content(f'<html><body style="margin:0">{svg}</body></html>')
        page.evaluate("document.fonts.ready")
        page.screenshot(path=str(output))
        browser.close()
    print(f"Upload {output} to the social_image_source configured in routes.toml")


if __name__ == "__main__":
    main()
