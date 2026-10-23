#!/bin/bash
#
# Build and release the BEAST 2.8 package.
#
# Usage:
#   ./release.sh              # run the tests and build the package ZIP in target/
#   ./release.sh --central    # also publish to Maven Central (mvn -Prelease deploy)
#   ./release.sh --github     # also create a GitHub release with the ZIP attached
#   ./release.sh --central --github
#
# Publishing requires a release version (no -SNAPSHOT) that agrees between pom.xml and
# version.xml, a "central" server with a Central Portal token in ~/.m2/settings.xml, a GPG
# signing key, and the gh CLI for --github.
#
set -euo pipefail
cd "$(dirname "$0")"

CENTRAL=false; GITHUB=false
for arg in "$@"; do
    case "$arg" in
        --central) CENTRAL=true ;;
        --github)  GITHUB=true ;;
        *) echo "Unknown option: $arg" >&2; exit 1 ;;
    esac
done

# --- Versions ---

PKG_LINE=$(grep '<package ' version.xml | head -1)
PKG_NAME=$(echo "$PKG_LINE" | sed "s/.*name=['\"]\\([^'\"]*\\)['\"].*/\\1/")
VERSION=$(echo "$PKG_LINE" | sed "s/.*version=['\"]\\([^'\"]*\\)['\"].*/\\1/")
POM_VERSION=$(mvn -q help:evaluate -Dexpression=project.version -DforceStdout)
POM_PKG_VERSION=$(mvn -q help:evaluate -Dexpression=beast.pkg.version -DforceStdout)

echo "=== ${PKG_NAME}: version.xml ${VERSION}, pom ${POM_VERSION} ==="
if [[ "$POM_PKG_VERSION" != "$VERSION" || "${POM_VERSION%-SNAPSHOT}" != "$VERSION" ]]; then
    echo "ERROR: version.xml (${VERSION}), pom version (${POM_VERSION}) and beast.pkg.version" \
         "(${POM_PKG_VERSION}) disagree" >&2
    exit 1
fi
if [[ ( "$CENTRAL" == true || "$GITHUB" == true ) && "$POM_VERSION" == *-SNAPSHOT ]]; then
    echo "ERROR: set the pom version to ${VERSION} (without -SNAPSHOT) before publishing" >&2
    exit 1
fi
if [[ ( "$CENTRAL" == true || "$GITHUB" == true ) && -n "$(git status --porcelain)" ]]; then
    echo "ERROR: commit all changes before publishing, so the release matches a commit" >&2
    exit 1
fi

# --- Build, test and assemble the package ZIP ---

mvn clean package
ZIP_NAME="${PKG_NAME}.package.v${VERSION}.zip"
ZIP_PATH="target/${ZIP_NAME}"
[[ -f "$ZIP_PATH" ]] || { echo "ERROR: expected ${ZIP_PATH}" >&2; exit 1; }
echo "=== Package: ${ZIP_PATH} ==="
unzip -l "$ZIP_PATH"

# --- Maven Central ---

if [[ "$CENTRAL" == true ]]; then
    echo "--- Publishing ${POM_VERSION} to Maven Central ---"
    mvn -Prelease deploy -DskipTests
    echo "Install with: packagemanager -maven io.github.ewanciuffi:multitypespikemodel:${VERSION}"
fi

# --- GitHub release ---

if [[ "$GITHUB" == true ]]; then
    GITHUB_REPO=$(git remote get-url origin \
        | sed 's|.*github.com[:/]\(.*\)\.git$|\1|; s|.*github.com[:/]\(.*\)$|\1|')
    echo "--- Creating GitHub release v${VERSION} on ${GITHUB_REPO} ---"
    gh release create "v${VERSION}" "$ZIP_PATH" \
        --repo "$GITHUB_REPO" \
        --target "$(git rev-parse HEAD)" \
        --title "${PKG_NAME} v${VERSION} (BEAST 2.8)" \
        --generate-notes
fi
