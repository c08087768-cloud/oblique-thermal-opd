# AGENTS.md — LaTeX Manuscript Workflow

## 1. Purpose and Roles

This repository contains an academic manuscript written in LaTeX.

The goal is to maintain a scientifically accurate, consistent, reproducible, and well-structured manuscript while keeping all changes traceable and easy to review.

### Responsibilities

- **The user** makes the final decisions about scientific content, manuscript structure, and unresolved academic questions.
- **ChatGPT** is used to discuss research logic, scientific explanations, mathematical derivations, and academic English. Confirmed content is handed to Codex as LaTeX-ready text with explicit editing instructions.
- **Codex** implements approved changes, checks manuscript-wide consistency, validates the LaTeX project, and manages version control according to this document.
- **GitHub** stores the revision history and supports synchronization, review, and recovery.

Do not treat a discussion, suggestion, or unconfirmed draft as approved manuscript content.

## 2. Editing Permissions

Use limited autonomy.

### Changes Codex may make autonomously

- Correct obvious spelling mistakes.
- Fix mechanical LaTeX formatting errors.
- Fix unambiguous broken paths or reference errors when the correct target can be established.
- Correct clearly identifiable, non-substantive caption or formatting mistakes.
- Perform necessary compilation and validation checks.

### Changes requiring prior confirmation

Report the issue and wait for approval before making changes that may affect:

- Scientific meaning or physical interpretation.
- Mathematical derivations, equations, or boundary conditions.
- Definitions of physical quantities or mathematical symbols.
- Experimental data, numerical results, or conclusions.
- Academic arguments, paragraph meaning, or section structure.
- Substantive English rewriting or stylistic polishing.

Do not silently resolve scientific inconsistencies.

If an inconsistency directly affects the correctness of the requested change, pause that part of the task and ask for clarification. Unrelated, safe work may continue.

## 3. ChatGPT-to-Codex Handoff

When implementing content discussed with ChatGPT:

1. Use the explicitly confirmed LaTeX-ready content and editing instructions provided for the task.
2. Locate the relevant source file and insertion or replacement point.
3. Inspect the existing macros, notation, equation labels, citation keys, and document conventions.
4. Preserve the approved scientific meaning and wording. Do not independently rewrite confirmed text.
5. Integrate the content into the existing manuscript structure.
6. Check related material elsewhere in the manuscript.
7. Compile and validate the changes when possible.
8. Commit and synchronize the work according to the Git workflow below.

Do not invent citation keys, equation labels, or macro definitions. If a required identifier or project convention cannot be established, report the uncertainty rather than guessing.

## 4. Git and GitHub Workflow

### Branch policy

- `writing` is the working branch for manuscript development and routine synchronization.
- `main` represents the stable, reviewed manuscript version.
- Do not merge into `main` without the user's explicit approval.

### Before modifying files

- Inspect the current branch and Git working-tree status.
- Identify existing uncommitted changes.
- Preserve the user's work; do not overwrite or discard unrelated changes.
- Inspect the repository's existing remotes and branch configuration before performing Git operations.

### During development

- Make small, logically related changes.
- After each meaningful, independently reviewable change, run the applicable checks, create a commit, and push to `writing`.
- Use descriptive commit messages.
- Do not create commits for changes that were not actually made.
- Never use destructive Git operations, force-push, or discard user changes without explicit authorization.

### When compilation fails

- Compilation failure does not prevent saving work to `writing`.
- Commit and push the work when appropriate, while explicitly recording the failure and unresolved issues.
- Do not claim that the manuscript has passed validation.
- Do not merge an unvalidated, failing version into `main`.

### Stable-version promotion

Before proposing a merge into `main`:

- Complete the required manuscript-wide checks.
- Resolve relevant compilation and reference errors.
- Report unresolved scientific or technical issues.
- Provide the commit identifier and validation status.
- Wait for the user's explicit approval before merging.

## 5. Manuscript Consistency Checks

Use three levels of review.

### Level 1 — Local checks

After every modification, inspect the changed content and its immediate context.

Check LaTeX syntax, mathematical notation, citations, labels, and nearby prose.

### Level 2 — Cross-section checks

When a change affects a symbol, definition, technical term, equation, assumption, or scientific claim, search the entire manuscript for related occurrences.

Include relevant:

- Equations and mathematical expressions.
- Main text and section headings.
- Figure and table captions.
- Cross-references and related explanations.
- LaTeX source files included by the main document.

Identify affected locations and report inconsistencies.

Do not perform a blind global replacement merely to make notation look uniform. Verify the meaning of each occurrence before changing it.

### Level 3 — Full manuscript review

Before proposing a merge into `main`, review the manuscript as a whole.

Check:

- Consistency of notation and terminology.
- Agreement between assumptions, derivations, results, and conclusions.
- Figure and table references.
- Equation labels and cross-references.
- Citation completeness and bibliography consistency.
- LaTeX compilation and other applicable validation results.

Report scientific inconsistencies instead of silently resolving them.

## 6. LaTeX Source Formatting

### Paragraphs

Use semantic line breaks in manuscript source files:

- Put each English sentence on a separate source line.
- Separate paragraphs with one blank line.
- Do not use `\\` to force ordinary prose line breaks.
- Avoid unrelated reformatting of unchanged paragraphs.

These rules concern source readability and version-control diffs, not the final PDF line wrapping.

### Equations

- Number only important equations or equations that need to be referenced.
- Let LaTeX generate equation numbers automatically.
- Never type equation numbers manually.
- Use appropriate numbered or unnumbered mathematical environments.
- Give referenced equations unique labels using the existing naming convention.
- Reference equations using `\eqref{...}` when compatible with the project.
- Use aligned mathematical environments for multistep derivations where appropriate.
- Do not change mathematical meaning merely to improve formatting.
- Do not remove existing equation labels solely for stylistic consistency.

Follow the existing document class, journal template, macros, and package configuration.

## 7. Citations and Bibliography

- Follow the citation and bibliography system already configured in the project.
- Never fabricate papers, authors, DOIs, publication years, or bibliographic metadata.
- Preserve existing citation keys unless a change is necessary and all affected references can be updated safely.
- Check for undefined citation keys and missing references after editing.
- Verify bibliographic metadata before adding new entries.
- Do not treat a candidate reference as approved evidence without the user's confirmation.

When a claim appears to require additional support:

1. Identify the specific claim.
2. Explain what kind of evidence is needed.
3. Search for and verify relevant candidate sources when appropriate.
4. Report the candidates and their relevance.
5. Wait for approval before incorporating them as accepted supporting evidence.

Do not add citations merely because a paper appears topically related.

## 8. Figures, Tables, and Experimental Results

- Check figure and table paths, labels, references, captions, and units.
- Check consistency between captions, surrounding prose, and the actual figures or tables.
- When a mismatch is suspected, inspect the available source figure, plotting scripts, data, and relevant context before deciding what is wrong.
- Correct an unambiguous caption, path, or reference mistake when the correction does not alter scientific meaning.
- Report unresolved discrepancies involving physical interpretation, experimental results, or conclusions.
- Do not modify raw experimental data, numerical results, or source images without explicit authorization.
- Do not regenerate or substantively alter scientific figures without explicit task authorization.

Never change data merely to make a figure agree with the manuscript text.

## 9. Compilation and Validation

- Identify and follow the existing build procedure.
- Use the compiler and build configuration already established by the project.
- Compile after changes whenever the environment permits.
- Check for compilation errors, missing figures, undefined references, duplicate labels, and relevant warnings.
- Distinguish newly introduced problems from pre-existing errors whenever possible.
- Do not delete manuscript content, suppress errors, or change scientific content merely to obtain a successful build.
- Report checks that were not performed.

Compilation success does not by itself establish scientific correctness.

If compilation fails, report the failure and preserve the work on `writing` when appropriate. Do not present an unvalidated version as stable.

## 10. Completion Report

After completing each task, provide a concise report containing the following five items.

### 1. Changes made

List the files and sections modified and summarize the work performed.

### 2. Scientific impact

State whether equations, notation definitions, physical assumptions, experimental data, numerical results, or conclusions were changed.

### 3. Validation results

Report compilation status, consistency checks, reference checks, and any checks that could not be performed.

### 4. GitHub synchronization

Report the working branch, commit identifier, push status, and any remaining uncommitted changes.

### 5. Outstanding issues

List unresolved errors, scientific questions requiring the user's decision, and necessary follow-up actions.

Do not claim success for operations that were not completed.

By default, provide this report in the Codex conversation. Do not create additional report files unless requested or necessary for a specific long-term record.

## 11. General Change Discipline

- Prefer small, targeted, reviewable changes.
- Preserve existing project structure unless a structural change is explicitly requested.
- Do not combine unrelated manuscript edits or code refactoring into one task.
- Do not silently expand the scope of an editing request.
- Report meaningful inconsistencies even when they are outside the immediate editing scope.
- Keep changes traceable through Git commits.
- When uncertain about scientific meaning or authorization, ask rather than guess.

The priority is to preserve scientific correctness while making manuscript development efficient, understandable, and reversible.