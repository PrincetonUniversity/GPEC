# Claude Coding Guidelines for GPEC

## Fortran Fixed-Format (.f files) Requirements

**CRITICAL: Line Width Limit**
- Fixed-format Fortran files (.f, NOT .f90) have a **72-character line limit**
- Lines exceeding 72 characters will cause compilation errors or be silently truncated
- **ALWAYS** check line length before making edits to .f files

### Line Format Rules
- Columns 1-5: Statement labels (optional)
- Column 6: Continuation character (use `$` or `&` for continuation lines)
- Columns 7-72: Fortran statements
- Columns 73-80: Ignored (historically used for sequence numbers)

### Continuation Lines
When a statement exceeds 72 characters, split it across multiple lines:
```fortran
      variable_name = very_long_expression_that_would_exceed_limit
     $                + more_of_the_expression
```

### String Splitting
For long strings, use Fortran string concatenation:
```fortran
      WRITE(*,*) "This is a very long warning message that "//
     $           "needs to be split across multiple lines"
```

## Type Precision

**REAL Variables:**
- Use `0.0_r8` instead of `0.0` for REAL(r8) variables
- Use `1.0_r8` instead of `1.0` for REAL(r8) variables
- Type suffix `_r8` ensures proper precision matching

## Variable Naming Conventions

Use consistent naming patterns within subroutines:
- `hw_*` prefix for half-widths (e.g., `hw_isl`, `hw_v`, `hw_v_crit`, `hw_sat`, `hw_min`)
- Avoid mixing naming styles within the same subroutine

## Safety Checks

**Array Access:**
- Always check if allocatable arrays are ALLOCATED before accessing them:
```fortran
      IF (ALLOCATED(array_name)) THEN
         ! Use array
      ELSE
         ! Handle unallocated case
      ENDIF
```

**Division:**
- Guard against division by zero for potentially small values
- Check bounds for ACOS, ASIN arguments (must be in [-1, 1])
- Check that square root arguments are non-negative

## Before Committing

**Pre-commit Checklist:**
1. Check all modified .f files for lines exceeding 72 characters
2. Verify type suffixes match variable declarations (_r8 for REAL(r8))
3. Ensure continuation lines use proper column 6 markers ($)
4. Test compilation before committing
5. Ask user permission before committing changes

## Verification Commands

Check for lines exceeding 72 characters:
```bash
awk 'length > 72 {print NR": "length" chars"}' filename.f
```

Count characters in a specific line:
```bash
sed -n 'NUMp' filename.f | wc -c
```
