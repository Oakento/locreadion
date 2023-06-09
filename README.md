# locreadion
This is to remove the ambiguity in counting the number of reads mapped to regions, as some reads could be counted multiple times if they span cross regions or include splice sites.

## Requirements
- bedtools
- chromosome name starts with `chr`

## Usage
```bash
locreadion -a <alignment> -r <region_dir> -o <output_dir>
```
