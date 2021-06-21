# ClassPro

A k-mer classification program using [FASTK](https://github.com/thegenemyers/FASTK) and [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) libraries by Gene Myers.

```
Classifier [-vs] [-T<int(4)>] [-c<int>] [-r<int>] <source_root>
```

- Input: `<source_root>.db`, `<source_root>.hist`, and `<source_root>.prof` made from the same input read dataset.
- Output: A fastq-like file `<source_root>.class` and a DAZZ_DB track `.<source_root>.class.[anno|data]`as a DAZZ_DB track (0: Error, 1: Repeat, 2: Haplo, 3: Diplo).
