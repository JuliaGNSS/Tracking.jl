# Agent instructions

## Commit messages: Conventional Commits (drive the releases)

Releases are fully automated: every push to `master` runs semantic-release
(`.github/workflows/Release.yml`), which derives the version bump, the
changelog, and the JuliaRegistrator call from the commit messages. Rebase
merge is enforced, so **every commit in a PR lands verbatim on `master`** and
must follow [Conventional Commits](https://www.conventionalcommits.org):

```
type(scope): short imperative summary

optional body explaining what and why
```

- `feat` → minor release, `fix` / `perf` → patch release; `refactor`, `test`,
  `docs`, `build`, `chore` → no release.
- Pick the scope from the area touched (e.g. `track`, `dc`, `onebit`,
  `bench`); look at `git log --oneline` for precedent.

### Breaking changes must be marked and explained

A breaking change **must** carry both:

1. a `!` after the type/scope — e.g. `feat(track)!: …`, and
2. a `BREAKING CHANGE:` footer paragraph that spells out *what* breaks and
   *how to migrate* (removed/renamed exported types, functions, fields, or
   kwargs; changed constructor signatures; changed observable behavior of
   public API).

The footer is what triggers the major version bump and becomes the release
notes — an unmarked breaking change silently ships as a minor/patch release.
When in doubt: it is breaking if released user code — including custom
`AbstractDopplerEstimator` implementations and direct
`downconvert_and_correlate(!)` callers — could stop compiling or change
behavior after the release.

### Do not hand-maintain release artifacts

The changelog and the version are produced by CI: on `master`,
semantic-release generates `CHANGELOG.md` entries and the release notes from
the commit messages and bumps `Project.toml`'s `version` accordingly. Never
edit `CHANGELOG.md` or bump `version` manually in a PR — write the commit
message well instead, since that text *is* the changelog.

## Formatting

CI enforces JuliaFormatter with the exact version pinned in
`.github/workflows/format.yml` (currently 2.8.5) against
`.JuliaFormatter.toml`. Before committing:

```julia
using JuliaFormatter; format(["src", "test"])
```
