# CONTRIBUTING

Contributions are always welcome.

## Formatting

- Julia code: format with [Runic](https://github.com/fredrikekre/Runic.jl).
  Keep lines under ~92 characters.
- Markdown: format with [Prettier](https://prettier.io):
  `npx prettier --write "**/*.md"`.
  Keep lines under ~100 characters.
  (URLs, LaTeX, code blocks, and tables may exceed these limits
  when longer lines read better.)

## Testing

Run the test suite with

```sh
julia --project=. -e 'using Pkg; Pkg.test("RAS_DMFT")'
```

## Commits

Follow [Conventional Commits](https://www.conventionalcommits.org/en/v1.0.0/);
see the git history for examples.
