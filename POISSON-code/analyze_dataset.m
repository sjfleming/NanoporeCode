function analyze_dataset(name,ref)

    fpath = ['C:\Minion\' name];
    analyze_fasta(fpath);
    align_fasta(fpath,ref);
    analyze_bam(fpath);

end