//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

class Utils {

    //
    // When running with -profile conda, warn if channels have not been set-up appropriately
    //
    public static void checkCondaChannels(log) {
        Yaml parser = new Yaml()
        def channels = []
        try {
            def config = parser.load("conda config --show channels".execute().text)
            channels = config.channels
        } catch(NullPointerException | IOException e) {
            log.warn "Could not verify conda channel configuration."
            return
        }

        // Check that all channels are present
        // This channel list is ordered by required channel priority.
        def required_channels_in_order = ['conda-forge', 'bioconda', 'defaults']
        def channels_missing = ((required_channels_in_order as Set) - (channels as Set)) as Boolean

        // Check that they are in the right order
        def channel_priority_violation = false
        def n = required_channels_in_order.size()
        for (int i = 0; i < n - 1; i++) {
            channel_priority_violation |= !(channels.indexOf(required_channels_in_order[i]) < channels.indexOf(required_channels_in_order[i+1]))
        }

        if (channels_missing | channel_priority_violation) {
            log.warn "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n" +
                "  There is a problem with your Conda configuration!\n\n" +
                "  You will need to set-up the conda-forge and bioconda channels correctly.\n" +
                "  Please refer to https://bioconda.github.io/\n" +
                "  The observed channel order is \n" +
                "  ${channels}\n" +
                "  but the following channel order is required:\n" +
                "  ${required_channels_in_order}\n" +
                "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
        }
    }

    public static Map flowcellLaneFromFastq(path) {
        def line
	path.withInputStream {
	    InputStream gzipStream = new java.util.zip.GZIPInputStream(it)
	    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
	    BufferedReader buffered = new BufferedReader(decoder)
	    line = buffered.readLine()
        }
	assert line.startsWith('@')
	line = line.substring(1)
	line = line.split("/")[0]
	line = line.split("\\s")[0]
	def fields = line.split(':')
	String fcid
	String lane
	if (fields.size() >= 7) {
            fcid = fields[2]
	    lane = fields[0]
	} else if (fields.size() == 5) {
	    fcid = fields[0]
	    lane = fields[1]
	} else {
            fcid = fields[0]
	    lane = ""
	}
	return ["fcid":fcid, "lane":lane]
    }
}
