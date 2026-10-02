from __future__ import annotations

from pathlib import Path

from tools.audit_public_tree import audit_tree, main


def test_clean_code_tree_passes(tmp_path: Path) -> None:
    (tmp_path / "module.py").write_text(
        '"""Small, generic source module."""\nVALUE = 3\n',
        encoding="utf-8",
    )
    (tmp_path / "README.md").write_text(
        "# Generic method\n\nNo study data are distributed.\n",
        encoding="utf-8",
    )

    assert audit_tree(tmp_path) == []


def test_prohibited_data_extension_is_reported(tmp_path: Path) -> None:
    path = tmp_path / "sample.fasta"
    path.write_text("placeholder", encoding="utf-8")

    violations = audit_tree(tmp_path)

    assert any(item.path == Path("sample.fasta") for item in violations)


def test_unrecognized_binary_file_type_is_reported(tmp_path: Path) -> None:
    path = tmp_path / "opaque.bin"
    path.write_bytes(bytes(range(20)))

    violations = audit_tree(tmp_path)

    assert any(
        item.path == Path("opaque.bin")
        and item.reason == "file type is not allowed in the code-only public tree"
        for item in violations
    )


def test_sensitive_content_patterns_are_reported(tmp_path: Path) -> None:
    site_name = "Con" + "sett"
    home_path = "/" + "Users" + "/researcher/private"
    project_identifier = "PS-" + "A1B2C3D4"
    protein_like_text = "A" * 60
    (tmp_path / "notes.txt").write_text(
        "\n".join([site_name, home_path, project_identifier, protein_like_text]),
        encoding="utf-8",
    )

    reasons = {item.reason for item in audit_tree(tmp_path)}

    assert "contains the study-site name" in reasons
    assert "contains an absolute user-home path" in reasons
    assert "contains a project-specific protein identifier" in reasons
    assert "contains a long protein-like character string" in reasons


def test_secret_bearing_file_names_and_oversized_files_are_reported(
    tmp_path: Path,
) -> None:
    (tmp_path / ".env").write_text("placeholder", encoding="utf-8")
    (tmp_path / "large.txt").write_text("x" * 20, encoding="utf-8")

    violations = audit_tree(tmp_path, maximum_file_size_bytes=10)
    paths = {item.path for item in violations}

    assert Path(".env") in paths
    assert Path("large.txt") in paths


def test_placeholder_environment_example_is_allowed(tmp_path: Path) -> None:
    (tmp_path / ".env.example").write_text(
        "SERVICE_TOKEN=replace_with_local_value\n",
        encoding="utf-8",
    )

    assert audit_tree(tmp_path) == []


def test_generated_cache_directories_are_ignored(tmp_path: Path) -> None:
    cache = tmp_path / "__pycache__"
    cache.mkdir()
    (cache / "sample.fasta").write_text("placeholder", encoding="utf-8")

    assert audit_tree(tmp_path) == []


def test_command_returns_nonzero_for_a_violation(tmp_path: Path) -> None:
    (tmp_path / "model.hmm").write_text("placeholder", encoding="utf-8")

    assert main([str(tmp_path)]) == 1
