"""
Command line interface for svZeroDTuner.
"""

from __future__ import annotations

import argparse
import sys

from .sv0d_tuner import run_optimization
from .sensitivity import run_sensitivity_analysis


def _add_optimize_parser(subparsers: argparse._SubParsersAction) -> None:
    optimize_parser = subparsers.add_parser(
        "optimize",
        help="Run parameter optimization from a tuning config",
    )
    optimize_parser.add_argument(
        "config",
        help="Path to tuning YAML config file",
    )
    optimize_parser.set_defaults(_handler=_handle_optimize)

    run_parser = subparsers.add_parser(
        "run",
        help="Alias for optimize",
    )
    run_parser.add_argument(
        "config",
        help="Path to tuning YAML config file",
    )
    run_parser.set_defaults(_handler=_handle_optimize)


def _add_sensitivity_parser(subparsers: argparse._SubParsersAction) -> None:
    sensitivity_parser = subparsers.add_parser(
        "sensitivity-analysis",
        help="Run sensitivity analysis from a config",
    )
    sensitivity_parser.add_argument(
        "config",
        help="Path to sensitivity YAML config file",
    )
    sensitivity_parser.set_defaults(_handler=_handle_sensitivity)

    alias_parser = subparsers.add_parser(
        "sensitivity",
        help="Alias for sensitivity-analysis",
    )
    alias_parser.add_argument(
        "config",
        help="Path to sensitivity YAML config file",
    )
    alias_parser.set_defaults(_handler=_handle_sensitivity)


def _handle_optimize(args: argparse.Namespace) -> int:
    result = run_optimization(args.config)
    success = bool(result.get("success", False))
    if not success:
        message = result.get("message", "Optimization failed")
        print(f"[svzerodtuner] {message}")
        return 1
    return 0


def _handle_sensitivity(args: argparse.Namespace) -> int:
    try:
        run_sensitivity_analysis(args.config)
    except Exception as exc:
        print(f"[svzerodtuner] Sensitivity analysis failed: {exc}")
        return 1
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="svzerodtuner",
        description="svZeroDTuner command line interface",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)
    _add_optimize_parser(subparsers)
    _add_sensitivity_parser(subparsers)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    handler = getattr(args, "_handler", None)
    if handler is None:
        parser.print_help()
        return 2
    return handler(args)


if __name__ == "__main__":
    raise SystemExit(main())
