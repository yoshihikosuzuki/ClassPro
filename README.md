# FASTK_tools
Some programs using [FASTK](https://github.com/thegenemyers/FASTK) and [DAZZ_DB](https://github.com/thegenemyers/DAZZ_DB) libraries by Gene Myers.

```
Classifier <source_root> <diplo_depth<int>>
```

- Both `<source_root>.db` and `<source_root>.prof` must exist, and these should be made from homopolymer compressed sequences. Current implementation requires that `<source_root>.db` is a database of homopolymer compressed reads, if one uses hoco reads.
- `diplo_depth` is the global avegarage depth per diploid. It will be automatically determined in the near future.
- Outputs `.<sequence_root>.class.anno` and `.<sequence_root>.class.data` as a DAZZ_DB track (0: Error, 1: Repeat, 2: Haplo, 3: Diplo).

```
EmerRate [-t<int(5)>] <source_root>[.prof]
```

- For each read, it just calculates the number of k-mers whose counts are less than the threshold given by `-t` and outputs to stdout.
