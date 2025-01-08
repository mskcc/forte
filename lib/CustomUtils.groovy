//
// This file holds several Groovy functions that could be useful for any Nextflow pipeline
//

import org.yaml.snakeyaml.Yaml

class CustomUtils {

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
            lane = fields[3]
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
