# CONTRIBUTING

Contributions are always welcome.

## Formatting

- Julia code: format with [Runic](https://github.com/fredrikekre/Runic.jl):
  `runic --check .` reports, `runic -i .` rewrites.
  Keep lines under ~92 characters.
- Markdown: format with [Prettier](https://prettier.io):
  `prettier --write "**/*.md"`.
  Keep lines under ~100 characters.
- URLs, LaTeX, code blocks, tables, and hardcoded numbers in tests
  may exceed both limits when longer lines read better.
- Break Markdown lines at logical locations,
  such as commas, brackets, or the end of a sentence,
  and start every sentence on a new line.
  Prettier does not check this,
  as `proseWrap` keeps its default `preserve`,
  which is also why it must never be set to `always`.

## Naming

- Struct fields stay ASCII (`e_v`, `t_v`).
- Function arguments may use Unicode (`ϵ_imp`, `ϵ_mf`, `Δ`).
- Greek epsilon is U+03F5 `ϵ` (`\epsilon`),
  never U+03B5 `ε` (`\varepsilon`).
- Use the ASCII hyphen `-` in comments and docstrings,
  never U+2212 `−` or U+2013 `–`.

## Error messages

- Wrap an error message that interpolates a value in
  [`lazy"..."`](https://docs.julialang.org/en/v1/base/strings/#Base.LazyString):

  ```julia
  function check_positive(x::Int)
      x >= 0 || throw(ArgumentError(lazy"message $(x)"))
      return x
  end
  ```

- A message without interpolation stays a plain string literal,
  where `lazy"..."` would only add noise.

## Testing

Run the test suite from the repository root with

```sh
julia --project=. --eval 'using Pkg; Pkg.test()'
```

- A new feature comes with its tests.
- A bug fix comes with a test that fails on the unfixed version.
- Prefer hardcoded expected values over comparing a function to itself.
- Set a tolerance by rounding the worst deviation you measure up,
  otherwise the test goes flaky on another machine.
  Measure it at different thread counts,
  as threaded sums change the rounding,
  on every Julia version of the CI matrix,
  and on the CPU of the CI runner
  (`OPENBLAS_CORETYPE=Haswell julia -C haswell`),
  as SIMD width changes it too.
  Round 1.2 times the worst deviation up to one significant digit
  (`2`, `30`, `8.0e3` times `eps()`),
  and use a few ulp where you measured exactly zero.

## Commits

- Follow [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/).
  The git history has examples.
- The scope is optional and names the area a change is restricted to,
  whatever the type,
  as in `docs(Poles)` or `test(Combinatorics)`:
  - `Poles`
  - `Combinatorics`

  Leave it out as soon as the change reaches beyond one of them,
  as most of the history does.

- Keep one logical change per commit.
- Write the header only.
  Add a body when the diff alone would mislead,
  such as a breaking change,
  a non-obvious physics or numerics reason,
  or a fix whose cause is not visible in the changed lines.

## CHANGELOG

- Only user-visible changes get an entry.
  A refactor, a test, or a comment does not.
- Open the pull request first.
  Then add the entries in a trailing `chore: update CHANGELOG` commit
  and push again.
  The PR number does not exist before that,
  and a commit cannot contain its own hash.
- Cite both, as
  `([#NNN](https://github.com/tuwien-cms/RAS_DMFT.jl/pull/NNN)) (short hash)`.
- One PR contributes as many entries as it has changes.
  Each cites the commit that made its statement true.
