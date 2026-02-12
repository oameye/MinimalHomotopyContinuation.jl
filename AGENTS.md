# Formatting
- When writing Julia code, use _4 whitespaces_ for indentation and try to keep
  the maximum line length under _92 characters_.
- When writing Markdown text, use _2 whitespaces_ for indentation and try to
  keep the maximum line length under _80 characters_.
  - Additionally, prioritize simple text style and limit unnecessary decorations
    (e.g. `**`) to only truly necessary locations. This is a style that should
    generally be aimed for, but pay particular attention when writing Markdown.
  - Headers should use sentence case (only the first word capitalized), not
    title case. For example:
    - Good: `## Conclusion and alternative approaches`
    - Bad: `## Conclusion And Alternative Approaches`
- When writing commit messages, follow the format "component: Brief summary" for
  the title. In the body of the commit message, provide a brief prose summary of
  the purpose of the changes made.
  Use backticks for code elements (function names, variables, file paths, etc.)
  to improve readability.
  Also, ensure that the maximum line length never exceeds 72 characters.
  When referencing external GitHub PRs or issues, use proper GitHub interlinking
  format (e.g., "owner/repo#123" for PRs/issues).
  Finally, if you write code yourself, include a "Written by Claude" footer at
  the end of the commit message (no emoji nonsense). However, when simply asked
  to write a commit message, there's no need to add that footer.
- For file names, use `-` (hyphen) as the word separator by default.
  However, if the file name corresponds directly to Julia code (e.g., a module
  name), use `_` (underscore) instead, since Julia identifiers cannot contain
  hyphens (unless we use `var"..."`). For example, test files like
  `test_completions.jl` define a module `module test_completions`,
  so they use underscores.

# Coding rules
- When writing functions, use the most restrictive signature type possible.
  This allows JET to easily catch unintended errors.
  Of course, when prototyping, it's perfectly fine to start with loose type
  declarations, but for the functions you ultimately commit, it's desirable to
  use type declarations as much as possible.
  Especially when AI agents suggest code, please make sure to clearly specify
  the argument types that functions expect.
  In situations where there's no particular need to make a function generic, or
  if you're unsure what to do, submit the function with the most restrictive
  signature type you can think of.

- For function calls with keyword arguments, use an explicit `;` for clarity.
  For example, code like this:
  ```julia
  ...
  Position(; line=i-1, character=m.match.offset-1)
  ...
  ```
  is preferred over:
  ```julia
  ...
  Position(line=i-1, character=m.match.offset-1)
  ...
  ```

- For AI agents: **ONLY INCLUDE COMMENTS WHERE TRULY NECESSARY**.
  When the function name or implementation clearly indicates its purpose or
  behavior, redundant comments are unnecessary.

- On the other hand, for general utilities that expected to be used in multiple
  places in the language server, it's fine to use docstrings to clarify their
  behavior. However, even in these cases, if the function name and behavior are
  self-explanatory, no special docstring is needed.

# Running test code
Please make sure to test new code when you wrote.

When working on a specific component (e.g., completions, diagnostics),
run the component-specific test instead of the full test suite:
```bash
julia --startup-file=no --project=test -e 'using Test; @testset "test_XXX" include("test/test_XXX.jl")'
```
Note:
- `--startup-file=no` avoids loading unnecessary startup utilities
- `--project=test` enables the test env for proper test execution

For even faster iteration on a specific `@testset`, use
[TestRunner.jl](#using-testrunnerjl):
```bash
testrunner --project=test test/test_XXX.jl "testset_name"
```

Running `Pkg.test()` takes about 8 minutes (as of December 2025), so avoid it unless:
- Changes affect multiple components
- The user explicitly requests the full test suite
- You're unsure which tests are relevant

## Using TestRunner.jl

Additionally, by using `@testset` as shown above, not only are tests hierarchized,
but through integration with [TestRunner.jl](https://github.com/aviatesk/TestRunner.jl),
you can also selectively execute specific `@testset`s, without executing the
entire test file or test suite.
If you're using this language server for development as well, you can run tests
from code lenses or code actions within test files. If you need to run them from
the command line, you can use commands like the following
(assuming the `testrunner` executable is installed):
```bash
testrunner --project=test --verbose test/test_completions "some_completion_func"
```
Note that TestRunner.jl is still experimental.
The most reliable way to run tests is still to execute test files standalone.

# Environment-related issues
For AI agents: **NEVER MODIFY [Project.toml](./Project.toml) OR  [test/Project.toml](./test/Project.toml) BY YOURSELF**.
If you encounter errors that seem to be environment-related when running tests,
in most cases this is due to working directory issues, so first `cd` to the root directory of this project
and re-run the tests. Never attempt to fix environment-related issues yourself.
If you cannot resolve the problem, inform the human engineer and ask for instructions.

# About modifications to code you've written
If you, as an AI agent, add or modify code, and the user appears to have made
further manual changes to that code after your response, please respect those
modifications as much as possible.
For example, if the user has deleted a function you wrote, do not reintroduce
that function in subsequent code generation.
If you believe that changes made by the user are potentially problematic,
please clearly explain your concerns and ask the user for clarification.

# Git operations
Only perform Git operations when the user explicitly requests them.
After completing a Git operation, do not perform additional operations based on
conversational context alone. Wait for explicit instructions.

When the user provides feedback or points out issues with a commit:
- Do NOT automatically amend the commit or create a fixup commit
- Explain what could be changed, then wait for explicit instruction