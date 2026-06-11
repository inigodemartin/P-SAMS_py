# P-SAMS_py

Versión en Python de **P-SAMS** (Plant Small RNA Maker Suite): herramienta para
diseñar **amiRNA** (artificial microRNA) y **syn-tasiRNA**, prediciendo además
posibles dianas y off-targets mediante **TargetFinder**.

## Requisitos

- Python 3.9+ (probado con 3.12)
- [Biopython](https://biopython.org/) (solo necesario para construir bases de datos)
- Perl 5 (necesario para ejecutar `TargetFinder/targetfinder.pl`)
- `samtools` disponible en el `PATH` (se usa `samtools faidx`)

Los binarios de `ssearch36` (Linux, macOS arm64/x86_64) ya están incluidos en
`TargetFinder/` y deberían tener permisos de ejecución. Si no es así:

```bash
chmod +x TargetFinder/ssearch36_* TargetFinder/targetfinder.pl
```

## Instalación

```bash
git clone https://github.com/inigodemartin/P-SAMS_py.git
cd P-SAMS_py
pip install biopython
```

No hace falta ningún paso de "build" adicional: `psams.py` y `create_db.py` se
ejecutan directamente con `python3`.

---

## 1. Construir la base de datos de una especie

Antes de poder analizar genes de una especie con `-a/--accessions` hace falta
generar su base de datos de k-mers (usada para predecir off-targets) y registrarla
en `psams.conf`.

### Ficheros de entrada

Necesitas los ficheros típicos de **Phytozome** para tu especie:

- Un FASTA de transcritos (mRNA), p. ej. `Nicotiana_benthamiana.transcript.fa`
- Un fichero de anotación `*.annotation_info.txt` (formato Phytozome, 16 columnas
  separadas por tabulador)

### Comando

```bash
python3 create_db.py \
    -f Nicotiana_benthamiana.transcript.fa \
    -d Nicotiana_benthamiana.annotation_info.txt \
    -s Nicotiana_benthamiana \
    -v v1.0 \
    -k 15
```

Argumentos:

| Flag | Descripción |
|------|-------------|
| `-f, --fasta`        | FASTA de transcritos (formato Phytozome) |
| `-d, --descriptions` | Fichero de anotación (`annotation_info.txt`) |
| `-s, --species`      | Nombre/código de la especie (se usará luego con `-s` en `psams.py`) |
| `-v, --version`      | Versión del transcriptoma (solo para nombrar los ficheros generados) |
| `-k, --ksize`        | Tamaño de k-mer para la base de datos de off-targets. **Debe ser 15** (coincide con la semilla `SEED` usada por `psams.py`) |
| `-l, --lowmem`       | Opcional. Construye la base de datos directamente, más lento pero usando mucha menos RAM |

Este paso puede tardar bastante (recorre todo el transcriptoma generando k-mers).

### Qué genera

- `db/<especie>_<version>_prueba.db`: base de datos SQLite con las tablas `kmers`
  y `annotation`.
- Ficheros intermedios formateados junto al FASTA de entrada
  (`<especie>.<version>.transcripts.fasta`, `<especie>.<version>.annotation.txt`)
  y su índice `.fai`.
- Una nueva entrada en `psams.conf` (en la raíz del proyecto) con el formato:

```ini
[Nicotiana_benthamiana]
mRNA=/ruta/absoluta/a/Nicotiana_benthamiana.transcript.fa
sql=/ruta/absoluta/al/proyecto/db/Nicotiana_benthamiana_v1.0_prueba.db
```

`psams.conf` no está versionado (está en `.gitignore`); cada instalación tiene el
suyo con las especies que haya construido.

---

## 2. Ejecutar un análisis

Una vez la especie está registrada en `psams.conf`, se puede analizar un gen por
su accession:

```bash
python3 psams.py -a Nbe01g01610.7 -s Nicotiana_benthamiana -o runs/Nbe01g01610 -u
```

### Flags principales

| Flag | Descripción |
|------|-------------|
| `-a, --accessions`   | Accession(es) del gen, separados por comas. Requiere `-s` |
| `-f, --fasta`        | Alternativa a `-a`: fichero FASTA con la(s) secuencia(s) a analizar |
| `-s, --species`      | Especie tal y como aparece en `psams.conf`. Requerido si se usa `-a`, o si se usa `-f` y se quiere predicción de off-targets |
| `-o, --output_path`  | Carpeta donde se crean los resultados (por defecto, el directorio actual) |
| `-c, --construct`    | `amiRNA` (por defecto) o `syntasiRNA` |
| `-t, --foldback`     | `eudicot` (por defecto) o `monocot` |
| `-n, --noofftarget`  | Desactiva la predicción de off-targets con TargetFinder |
| `-u, --unlimit`      | No limitar a 3 resultados óptimos: recorre todos los candidatos posibles (más lento) |

Para `syntasiRNA`, se pueden definir varios grupos de genes/secuencias separando
los grupos con `:`, p. ej. `-a gen1,gen2:gen3,gen4`.

### Resultados generados

```
runs/Nbe01g01610/Nbe01g01610.7_psams_output/
├── Nbe01g01610.7_optimal_results.tsv      # amiRNAs óptimos (sin off-targets)
├── Nbe01g01610.7_suboptimal_results.tsv   # amiRNAs subóptimos (con off-targets)
├── Nbe01g01610.7_psams.json               # resultado final combinado
└── tf_results/
    ├── site_0001_TargetFinder_result.json # salida de TargetFinder por candidato
    └── ...
```

Si `-o` ya contiene resultados completos para esos accessions, el script lo detecta
y termina sin volver a ejecutar el análisis.

---

## 3. Prueba rápida (sin base de datos)

Para comprobar que la instalación funciona sin necesidad de construir ninguna base
de datos, puedes pasar tu propia secuencia con `-f` y desactivar la predicción de
off-targets con `-n` (este modo no usa Perl, samtools ni `psams.conf`).

Crea un fichero `example.fasta`:

```fasta
>my_transcript
ATGGCGGATTCAGAGAAGCCGGTTACCGGAAGCTTGAGCTCGGATCCACTAGTAACGGCCGCCAGTGTG
```

Y ejecútalo:

```bash
python3 psams.py -f example.fasta -n -o runs/quick_test
```

Esto generará `runs/quick_test/my_transcript_psams_output/` con hasta 3 diseños
de amiRNA óptimos (guide, star, oligos para clonación) calculados directamente
sobre la secuencia de entrada, sin comprobar off-targets.
