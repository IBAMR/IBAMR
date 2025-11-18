# Code Formatting Guide for IBAMR-Style C++

This guide explains how to maintain consistent code formatting that matches IBAMR coding standards.

## Quick Start

### Option 1: Using clang-format (Recommended)

If you have `clang-format` installed:

```bash
# Format a single file
clang-format -i -style=file IBNACA0012Kinematics.cpp

# Format all C++ files
clang-format -i -style=file *.cpp *.h *.C

# Preview changes without modifying files
clang-format -style=file IBNACA0012Kinematics.cpp | diff IBNACA0012Kinematics.cpp -
```

### Option 2: Manual Style Compliance

If `clang-format` is not available, follow the manual checklist below.

## IBAMR Style Requirements

### 1. Indentation
- **4 spaces** per indentation level
- **Never use tabs**
- Continuation lines: 4 spaces

```cpp
void IBNACA0012Kinematics::setKinematicsVelocity(const double time,
                                                  const std::vector<double>& incremented_angle,
                                                  const std::vector<double>& center_of_mass)
{
    // Function body indented 4 spaces
    if (condition)
    {
        // Nested block indented 4 more spaces
        doSomething();
    }
}
```

### 2. Braces
- **Opening brace** on the same line as function/control statement
- **Closing brace** on its own line, aligned with start of statement

```cpp
void function()
{
    // Good
}

class MyClass
{
public:
    // Good
};
```

### 3. Naming Conventions

| Element | Convention | Example |
|---------|-----------|---------|
| Classes | CamelCase | `IBNACA0012Kinematics` |
| Functions | camelCase | `setKinematicsVelocity()` |
| Member variables | `d_` + snake_case | `d_current_time` |
| Local variables | snake_case | `loop_time`, `iteration_num` |
| Constants | UPPER_CASE | `CHORD_LENGTH` |
| Namespaces | lowercase | `namespace ibamr` |

### 4. Spacing

**Around operators:**
```cpp
int x = 5 + 3;        // Good: spaces around binary operators
int y = func(a, b);   // Good: no space after ( or before )
```

**After keywords:**
```cpp
if (condition)        // Good: space after 'if'
for (int i = 0; i < n; ++i)  // Good: space after 'for'
```

**In function declarations:**
```cpp
void setShape(const double time, const std::vector<double>& angle);
// Good: no space before (, space after commas
```

### 5. Line Length
- **Preferred maximum:** 120 characters
- Break long lines at logical points

```cpp
// Good: Broken at natural point
Pointer<IBNACA0012Kinematics> ib_kinematics_op =
    new IBNACA0012Kinematics("naca0012carangiform",
                             app_initializer->getComponentDatabase("ConstraintIBKinematics"),
                             ib_method_ops->getLDataManager(),
                             patch_hierarchy);
```

### 6. Comments

**Doxygen style for public API:**
```cpp
/*!
 * \brief Set kinematics velocity for NACA0012 carangiform fish.
 *
 * \param time Current simulation time
 * \param incremented_angle_from_reference_axis Rotation angles
 * \param center_of_mass Center of mass position
 * \param tagged_pt_position Tagged point position
 */
void setKinematicsVelocity(const double time,
                           const std::vector<double>& incremented_angle_from_reference_axis,
                           const std::vector<double>& center_of_mass,
                           const std::vector<double>& tagged_pt_position);
```

**Regular comments:**
```cpp
// Single-line comments use // with space after

/*
 * Multi-line comments use this format
 * with asterisks aligned
 */
```

### 7. Include Order

```cpp
// 1. Config files
#include <SAMRAI_config.h>

// 2. External library headers
#include <petscsys.h>

// 3. SAMRAI headers
#include <BergerRigoutsos.h>
#include <CartesianGridGeometry.h>

// 4. IBAMR/IBTK headers
#include <ibamr/ConstraintIBMethod.h>
#include <ibtk/AppInitializer.h>

// 5. Application namespaces
#include <ibamr/app_namespaces.h>

// 6. Local headers
#include "IBNACA0012Kinematics.h"
```

### 8. Access Specifiers
- Not indented within class
- Separated by blank lines

```cpp
class MyClass
{
public:
    MyClass();
    void publicMethod();

private:
    void privateMethod();
    double d_member_variable;
};
```

### 9. Const Correctness
- Use `const` wherever possible
- Prefer `const&` for parameters

```cpp
// Good
void process(const std::string& name, const std::vector<double>& data);

// Member functions that don't modify state
const std::vector<std::vector<double>>& getShape(const int level) const;
```

### 10. Namespace Closing
- Add comment indicating which namespace

```cpp
namespace IBAMR
{
class MyClass
{
    // ...
};

} // namespace IBAMR
```

## Installation: clang-format

### Ubuntu/Debian
```bash
sudo apt-get install clang-format
```

### macOS
```bash
brew install clang-format
```

### Check version
```bash
clang-format --version
```

**Note:** IBAMR has historically used clang-format versions 6-8. Different versions may produce slightly different formatting.

## Pre-commit Hook (Optional)

To automatically format code before each commit:

```bash
# Create .git/hooks/pre-commit
cat > .git/hooks/pre-commit << 'EOF'
#!/bin/bash
# Format all staged C++ files
FILES=$(git diff --cached --name-only --diff-filter=ACM | grep -E '\.(cpp|h|C)$')
if [ -n "$FILES" ]; then
    clang-format -i -style=file $FILES
    git add $FILES
fi
EOF

chmod +x .git/hooks/pre-commit
```

## Verification Checklist

Before committing code, verify:

- [ ] Code compiles without warnings
- [ ] All files have IBAMR copyright header
- [ ] Include guards use `included_ClassName` pattern
- [ ] Proper include order (config → external → IBAMR → local)
- [ ] Classes use CamelCase, functions use camelCase
- [ ] Member variables have `d_` prefix
- [ ] 4-space indentation throughout
- [ ] No trailing whitespace
- [ ] Line length < 120 characters (preferred)
- [ ] Doxygen comments for public API
- [ ] Const-correctness applied
- [ ] Namespace closing comments present

## Reference

For the definitive IBAMR style:
- IBAMR Repository: https://github.com/IBAMR/IBAMR
- Official `.clang-format`: https://github.com/IBAMR/IBAMR/blob/master/.clang-format
- Formatting script: `IBAMR/scripts/reformat.sh`

## Current Status

Your code currently **passes all IBAMR style checks**. See `IBAMR_STYLE_COMPLIANCE_REPORT.md` for detailed analysis.

To maintain this standard:
1. Use the provided `.clang-format` file
2. Run `clang-format` before commits
3. Follow naming conventions consistently
4. Keep documentation up to date
