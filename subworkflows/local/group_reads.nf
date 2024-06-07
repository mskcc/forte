workflow GROUP_READS {

    take:
    ungrouped_reads

    main:

    grouped_reads = ungrouped_reads
        .map{ meta, reads ->
            def read_group = meta.read_group
            def fastq_pair_id = meta.fastq_pair_id
            def meta_clone = meta.clone().findAll { !["read_group","fastq_pair_id"].contains(it.key) }
            meta_clone.id = meta.sample
            [groupKey(meta_clone,meta.fq_num), reads, read_group, fastq_pair_id]
        }.groupTuple(by:[0])
        .map{ meta, reads, read_group, fastq_pair_id ->
            meta = meta + [read_group:read_group.join(','), fastq_pair_id:fastq_pair_id.join(',')]
            [meta, reads.flatten()]
        }

    emit:
    grouped_reads

}
