from typing import Any

from modflow_devtools.misc import try_get_enum_value


class Filters:
    @staticmethod
    def value(v: Any) -> str:
        """
        Format a value to appear in the RHS of an assignment or an argument-passing
        expression. Booleans become .true./.false.; everything else gets wrapped in
        single quotes.  Multi-line strings (produced by _wrap_f90_content in
        dfn2f90.py) are emitted as adjacent quoted strings joined with ``// &``
        Fortran concatenation so that each source line fits within the project
        line-length limit.
        """
        v = try_get_enum_value(v)
        if isinstance(v, bool):
            return ".true." if v else ".false."
        s = str(v)
        if "\n" in s:
            parts = s.split("\n")
            quoted = [f"'{p}'" for p in parts]
            return "// &\n    ".join(quoted)
        return f"'{v}'"
