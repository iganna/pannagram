#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import os
import re
import sys
from pathlib import Path

LOCUS_RE = re.compile(
    r"^locus_(\d+)(?:_aligned)?\.fasta$",
    re.IGNORECASE
)

def parse_input_genomes(input_txt: Path) -> list[str]:
    """
    input.txt содержит пути. Имя генома = basename пути.
    Порядок сохраняется. Дубликаты убираются (берём первый).
    """
    genomes = []
    seen = set()
    with input_txt.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            base = os.path.basename(line.rstrip("/\\"))
            if base and base not in seen:
                genomes.append(base)
                seen.add(base)
    return genomes

def index_locus_files(folder: Path) -> dict[int, Path]:
    """
    Возвращает mapping: locus_number -> filepath
    """
    idx: dict[int, Path] = {}
    for p in folder.iterdir():
        if not p.is_file():
            continue
        m = LOCUS_RE.match(p.name)
        if m:
            idx[int(m.group(1))] = p
    return idx

def read_one_line_fasta(path: Path) -> dict[str, str]:
    """
    FASTA, где после каждого заголовка >name ровно ОДНА строка последовательности.
    Берём только первую строку после заголовка (если вдруг встретятся лишние — игнорируем).
    """
    data: dict[str, str] = {}
    current = None
    expecting_seq = False

    with path.open("r", encoding="utf-8", errors="replace") as f:
        for raw in f:
            line = raw.rstrip("\r\n")
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].strip() or None
                expecting_seq = True
                continue
            if expecting_seq and current is not None:
                data[current] = line
                expecting_seq = False

    return data

def alignment_length(locus_map: dict[str, str]) -> int:
    for seq in locus_map.values():
        return len(seq)
    return 0

def safe_filename(name: str) -> str:
    return re.sub(r'[<>:"/\\|?*\x00-\x1F]', "_", name)

def main():
    ap = argparse.ArgumentParser(
        description=(
            "Переворот locus_*.fasta в файлы по геномам; "
            "пропущенные локусы -> пустая строка, "
            "пропущенные геномы в локусе -> строка гэпов нужной длины. "
            "Чтение локусов батчами."
        )
    )
    ap.add_argument("-i", "--input", default="input.txt", help="Файл со списком путей (basename = имя генома).")
    ap.add_argument("-d", "--dir", default=".", help="Папка, где лежат locus_*.fasta.")
    ap.add_argument("-o", "--out", default="out", help="Папка для выходных файлов по геномам.")
    ap.add_argument("-b", "--batch", type=int, default=1000, help="Размер батча локусов (по умолчанию 1000).")
    args = ap.parse_args()

    folder = Path(args.dir).resolve()
    input_txt = Path(args.input).resolve()
    out_dir = Path(args.out).resolve()

    if not input_txt.exists():
        print(f"ERROR: input.txt не найден: {input_txt}", file=sys.stderr)
        sys.exit(1)

    genomes = parse_input_genomes(input_txt)
    if not genomes:
        print("ERROR: input.txt не содержит валидных путей/имён.", file=sys.stderr)
        sys.exit(1)

    locus_idx = index_locus_files(folder)
    if not locus_idx:
        print(f"ERROR: в папке {folder} не найдено файлов вида locus_*.fasta", file=sys.stderr)
        sys.exit(1)

    max_locus = max(locus_idx.keys())
    out_dir.mkdir(parents=True, exist_ok=True)

    # Подготовим выходные файлы: очистим/создадим, чтобы дальше только дописывать.
    out_paths: dict[str, Path] = {}
    for genome in genomes:
        out_path = out_dir / f"{safe_filename(genome)}.txt"
        out_paths[genome] = out_path
        with out_path.open("w", encoding="utf-8", newline="\n") as out:
            pass

    batch_size = max(1, args.batch)
    missing_files_total = 0

    # Идём батчами по локусам
    for start in range(1, max_locus + 1, batch_size):
        end = min(max_locus, start + batch_size - 1)

        # 1) читаем батч локусов в память
        locus_maps: dict[int, dict[str, str] | None] = {}
        locus_lens: dict[int, int] = {}

        for i in range(start, end + 1):
            fp = locus_idx.get(i)
            if fp is None:
                locus_maps[i] = None
                locus_lens[i] = 0
                missing_files_total += 1
            else:
                mp = read_one_line_fasta(fp)
                locus_maps[i] = mp
                locus_lens[i] = alignment_length(mp)

        # 2) дописываем строки для каждого генома за этот батч
        for genome in genomes:
            out_path = out_paths[genome]
            with out_path.open("a", encoding="utf-8", newline="\n") as out:
                for i in range(start, end + 1):
                    mp = locus_maps[i]
                    if mp is None:
                        out.write("\n")  # пустая строка
                    else:
                        if genome in mp:
                            out.write(mp[genome] + "\n")
                        else:
                            L = locus_lens[i]
                            out.write(("-" * L if L > 0 else "") + "\n")

        # (необязательно) прогресс
        print(f"OK: обработан батч локусов {start}..{end}")

    print(f"OK: создано файлов: {len(genomes)}")
    print(f"Диапазон локусов: 1..{max_locus}")
    if missing_files_total:
        print(f"Пропущенных locus_*.fasta в диапазоне: {missing_files_total} (в этих позициях пустые строки у всех)")
    print(f"Выходная папка: {out_dir}")

if __name__ == "__main__":
    main()
