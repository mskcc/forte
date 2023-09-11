workflow BAIT_INPUTS {

    main:

    Channel.from(params.baits ? params.baits : [])
        .map{ba ->
            [
                ba.keySet().collect(),
                ba.keySet().collect{ba[it]["baits"]},
                ba.keySet().collect{ba[it]["targets"]}
            ]
        }
        .transpose()
        .set{baits}

    emit:
    baits = baits

}
