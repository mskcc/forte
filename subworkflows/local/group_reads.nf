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
            [meta_clone, reads, read_group, fastq_pair_id]
        }.groupTuple(by:[0])
        .map{ meta, reads, read_group, fastq_pair_id ->
            // sort all the lists, but all with respect to the fastq_pair_id
            def zipped = fastq_pair_id.indices.collect { i ->
                [fastq_pair_id[i],reads[i],read_group[i]]
            }
            zipped.sort { a, b -> a[0] <=> b[0] }
            def sorted_fastq_pair_id = zipped.collect { it[0] }
            def sorted_reads         = zipped.collect { it[1] }
            def sorted_read_group    = zipped.collect { it[2] }
            def meta_clone = [:]
            meta_clone = meta + [read_group:sorted_read_group.join(','), fastq_pair_id:sorted_fastq_pair_id.join(',')]
            [meta_clone, sorted_reads.flatten()]
        }

    emit:
    grouped_reads

}
