from typing import Any

from modflow_devtools.misc import try_get_enum_value


class Filters:
    @staticmethod
    def value(v: Any) -> str:
        """
        Format a value to appear in the RHS of an assignment or an argument-passing
        expression. Booleans become .true./.false.; everything else gets wrapped in
        single quotes.
        """
        v = try_get_enum_value(v)
        if isinstance(v, bool):
            return ".true." if v else ".false."
        return f"'{v}'"
