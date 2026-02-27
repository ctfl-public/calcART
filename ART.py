"""
Export everything from the real implementation modules, 
regardless of how the user’s Python path is set up.
"""

try:
    # Preferred when importing from parent dir as a package.
    from calcART.calcART import *  # type: ignore[F403]
    from calcART.objects import *  # type: ignore[F403]
except ImportError:
    # Fallback when repo root itself is on sys.path.
    from calcART import *  # type: ignore[F403]
    from objects import *  # type: ignore[F403]
