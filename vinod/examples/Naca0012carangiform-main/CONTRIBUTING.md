# Contributing to NACA0012 Carangiform

Thank you for your interest in contributing to this project! This document provides guidelines for contributing to the NACA0012 undulatory swimming simulation project.

## Code of Conduct

- Be respectful and constructive in discussions
- Focus on what is best for the community
- Show empathy towards other contributors

## How to Contribute

### Reporting Bugs

Before creating a bug report, please check existing issues to avoid duplicates.

When filing a bug report, include:
- A clear, descriptive title
- Steps to reproduce the issue
- Expected vs. actual behavior
- Your environment (OS, IBAMR version, compiler version)
- Relevant log output or error messages
- If applicable, input files and simulation parameters

### Suggesting Enhancements

Enhancement suggestions are welcome! Please include:
- A clear description of the proposed feature
- The motivation and use case
- Example usage or pseudocode if applicable
- Any potential drawbacks or alternatives considered

### Pull Requests

1. **Fork the repository** and create your branch from `main`
2. **Follow the coding style** (see IBAMR_STYLE_COMPLIANCE_REPORT.md)
3. **Test your changes** thoroughly
4. **Update documentation** if needed
5. **Write clear commit messages** following this format:
   ```
   Short (50 chars or less) summary

   More detailed explanatory text, if necessary. Wrap at 72 characters.
   Explain the problem this commit solves and why you chose this approach.
   ```

## Development Workflow

### Setting Up Development Environment

1. Install IBAMR and dependencies (see INSTALL.md)
2. Clone your fork:
   ```bash
   git clone https://github.com/YOUR_USERNAME/Naca0012carangiform.git
   cd Naca0012carangiform
   ```
3. Add upstream remote:
   ```bash
   git remote add upstream https://github.com/vinodthale/Naca0012carangiform.git
   ```

### Making Changes

1. Create a feature branch:
   ```bash
   git checkout -b feature/your-feature-name
   ```
2. Make your changes
3. Format code according to IBAMR style:
   ```bash
   clang-format -i *.cpp *.h
   ```
4. Test your changes:
   ```bash
   mkdir build && cd build
   cmake ..
   make
   # Run simulations to verify
   ```
5. Commit and push:
   ```bash
   git add .
   git commit -m "Your descriptive commit message"
   git push origin feature/your-feature-name
   ```

### Coding Standards

This project follows IBAMR coding standards:

- **Indentation**: 4 spaces (no tabs)
- **Line length**: Maximum 120 characters
- **Naming conventions**:
  - Classes: `PascalCase` (e.g., `IBNACA0012Kinematics`)
  - Functions: `camelCase` (e.g., `getShape`)
  - Variables: `snake_case` (e.g., `time_step`)
  - Constants: `UPPER_SNAKE_CASE` (e.g., `MAX_ITERATIONS`)
- **Comments**: Use Doxygen-style comments for functions and classes
- **Headers**: Use `#pragma once` or include guards

See `.clang-format` for detailed formatting rules.

### Testing

Before submitting a pull request:

1. Verify mesh generation works for both modes:
   ```matlab
   swimming_mode = 'carangiform'; run('naca0012_swimmer_generator.m')
   swimming_mode = 'anguilliform'; run('naca0012_swimmer_generator.m')
   ```

2. Test simulation runs without errors:
   ```bash
   mpirun -np 4 ./main2d input2d
   ```

3. Check for memory leaks (if modifying C++ code):
   ```bash
   valgrind --leak-check=full ./main2d input2d
   ```

### Documentation

When adding new features:

- Update README.md with usage examples
- Add comments to complex algorithms
- Update MESH_PARAMETER_GUIDE.md if changing mesh generation
- Include references to relevant papers

## Areas for Contribution

We particularly welcome contributions in:

- **New swimming modes**: Subcarangiform, thunniform, ostraciiform
- **Performance optimization**: Parallel efficiency, memory usage
- **Validation**: Comparison with experimental data
- **Visualization**: Post-processing scripts for analysis
- **Testing**: Unit tests, regression tests
- **Documentation**: Tutorials, examples, API documentation
- **Bug fixes**: Issues tagged as "good first issue"

## Questions?

If you have questions about contributing:
- Open an issue with the "question" label
- Check existing documentation (README.md, INSTALL.md)
- Review IBAMR documentation: https://github.com/IBAMR/IBAMR

## Attribution

Contributors will be acknowledged in:
- Git commit history
- Release notes
- CONTRIBUTORS file (if significant contributions)

Thank you for contributing!
