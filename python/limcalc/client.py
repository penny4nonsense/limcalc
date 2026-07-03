"""
Client for the limcalc Haskell CLI.

Manages a persistent subprocess running limcalc-cli, communicates via
JSON line protocol (one JSON object per line in, one per line out),
and exposes a clean Python API for symbolic calculus operations.

The CLI binary is located by searching:
1. The LIMCALC_CLI environment variable
2. A bundled binary in the package directory
3. The system PATH

Example:
    from limcalc.client import LimCalc

    lc = LimCalc()
    result = lc.diff({"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}, "x")
    # {"tag": "Cos", "arg": {"tag": "Var", "name": "x"}}
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any


def _find_cli() -> str:
    """Locate the limcalc-cli executable.

    Search order:
    1. LIMCALC_CLI environment variable
    2. Bundled binary alongside this file (for PyPI wheels)
    3. System PATH

    Raises:
        FileNotFoundError: if no binary is found.
    """
    # 1. Environment variable override
    env_path = os.environ.get("LIMCALC_CLI")
    if env_path and Path(env_path).is_file():
        return env_path

    # 2. Bundled binary next to this file
    here = Path(__file__).parent
    for name in ["limcalc-cli", "limcalc-cli.exe"]:
        bundled = here / "bin" / name
        if bundled.is_file():
            return str(bundled)

    # 3. System PATH
    cli = shutil.which("limcalc-cli")
    if cli:
        return cli

    raise FileNotFoundError(
        "Could not find limcalc-cli executable. "
        "Set the LIMCALC_CLI environment variable to its path, "
        "or ensure limcalc-cli is on your PATH."
    )


class LimCalcError(Exception):
    """Raised when the limcalc CLI returns an error response."""
    pass


class LimCalc:
    """Client for the limcalc Haskell CLI.

    Manages a persistent subprocess and exposes symbolic calculus
    operations as Python methods. Each method sends a JSON command
    to the CLI and returns the parsed JSON result.

    Args:
        cli_path: Path to the limcalc-cli executable. If None,
                  uses _find_cli() to locate it automatically.

    Example:
        lc = LimCalc()
        result = lc.diff(
            {"tag": "Sin", "arg": {"tag": "Var", "name": "x"}},
            "x"
        )
    """

    def __init__(self, cli_path: str | None = None):
        self._cli_path = cli_path or _find_cli()
        self._proc: subprocess.Popen | None = None

    def _ensure_running(self) -> None:
        """Start the CLI subprocess if not already running."""
        if self._proc is None or self._proc.poll() is not None:
            self._proc = subprocess.Popen(
                [self._cli_path],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                encoding="utf-8",
            )

    def _send(self, command: dict) -> Any:
        """Send a command to the CLI and return the result.

        Args:
            command: A dict representing a limcalc CLI command.

        Returns:
            The parsed JSON result value.

        Raises:
            LimCalcError: if the CLI returns an error response.
            RuntimeError: if the CLI process dies unexpectedly.
        """
        self._ensure_running()
        line = json.dumps(command) + "\n"
        self._proc.stdin.write(line)
        self._proc.stdin.flush()

        response_line = self._proc.stdout.readline()
        if not response_line:
            stderr = self._proc.stderr.read()
            raise RuntimeError(
                f"limcalc-cli process died unexpectedly.\n"
                f"stderr: {stderr}"
            )

        response = json.loads(response_line)
        if not response.get("ok"):
            raise LimCalcError(response.get("error", "Unknown error"))

        return response["result"]

    def diff(self, expr: dict, var: str) -> dict:
        """Symbolically differentiate expr with respect to var.

        Args:
            expr: An Expr JSON object.
            var:  The variable name to differentiate with respect to.

        Returns:
            An Expr JSON object representing the derivative.

        Example:
            lc.diff({"tag": "Sin", "arg": {"tag": "Var", "name": "x"}}, "x")
            # {"tag": "Cos", "arg": {"tag": "Var", "name": "x"}}
        """
        return self._send({"op": "diff", "expr": expr, "var": var})

    def partial_diff(self, expr: dict, var: str) -> dict:
        """Partial derivative of expr with respect to var.

        Uses the algebraic chain-rule path (deriveExpr), which is
        correct for iterated partial derivatives.

        Args:
            expr: An Expr JSON object.
            var:  The variable name to differentiate with respect to.

        Returns:
            An Expr JSON object.
        """
        return self._send({"op": "partial_diff", "expr": expr, "var": var})

    def integrate(self, expr: dict, var: str) -> dict:
        """Symbolically integrate expr with respect to var (Risch algorithm).

        Args:
            expr: An Expr JSON object.
            var:  The variable name to integrate with respect to.

        Returns:
            An Expr JSON object representing the antiderivative.

        Raises:
            LimCalcError: if the integral is non-elementary or not implemented.
        """
        return self._send({"op": "integrate", "expr": expr, "var": var})

    def limit(self, expr: dict, var: str, x0: float) -> float:
        """Compute lim_{var -> x0} expr.

        Args:
            expr: An Expr JSON object.
            var:  The variable name.
            x0:   The limit point.

        Returns:
            The limit value as a float.

        Raises:
            LimCalcError: if the limit is a pole, does not exist, or
                          the expansion fails.
        """
        return self._send({"op": "limit", "expr": expr, "var": var, "x0": x0})

    def simplify(self, expr: dict) -> dict:
        """Apply algebraic simplification to expr.

        Args:
            expr: An Expr JSON object.

        Returns:
            A simplified Expr JSON object.
        """
        return self._send({"op": "simplify", "expr": expr})

    def pretty(self, expr: dict) -> str:
        """Pretty-print an Expr as a human-readable string.

        Args:
            expr: An Expr JSON object.

        Returns:
            A string like "2*x*cos(x^2)".
        """
        return self._send({"op": "pretty", "expr": expr})

    def gradient(self, expr: dict, vars: list[str]) -> list[dict]:
        """Compute the gradient of expr with respect to vars.

        Args:
            expr: An Expr JSON object (scalar-valued).
            vars: List of variable names.

        Returns:
            A list of Expr JSON objects [df/dx1, df/dx2, ...].
        """
        return self._send({"op": "gradient", "expr": expr, "vars": vars})

    def close(self) -> None:
        """Shut down the CLI subprocess."""
        if self._proc and self._proc.poll() is None:
            self._proc.stdin.close()
            self._proc.wait()
            self._proc = None

    def __enter__(self) -> "LimCalc":
        return self

    def __exit__(self, *args) -> None:
        self.close()
