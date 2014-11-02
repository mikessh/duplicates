/**
 Copyright 2014 Mikhail Shugay (mikhail.shugay@gmail.com)

 Licensed under the Apache License, Version 2.0 (the "License");
 you may not use this file except in compliance with the License.
 You may obtain a copy of the License at

 http://www.apache.org/licenses/LICENSE-2.0

 Unless required by applicable law or agreed to in writing, software
 distributed under the License is distributed on an "AS IS" BASIS,
 WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 See the License for the specific language governing permissions and
 limitations under the License.
 */

package com.antigenomics.duplicates

@Grab(group = 'org.codehaus.gpars', module = 'gpars', version = '1.2.1')

import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger
import java.util.concurrent.atomic.AtomicIntegerArray
import java.util.zip.GZIPInputStream

/**
 * Input arguments:
 *
 * .FASTQ[.gz] file name list
 * Those should initially come from migec/Checkout de-multiplexing utility,
 * should be post-processed using migec/CdrBlast with --cdr3-fastq-file,
 * having UMI:NNN:QQQ and CDR3:TGT... in their header
 * where NNN is the UMI sequence, QQQ is its quality Phred score string and TGT... is the sequence of CDR3
 *
 * Output file prefix (coming last)
 */

def inputFileNames = args[0..-2], outputPrefix = args[-1]

// Core data structure

class IlmnInfo {
    private final int tile, x, y, qSum
    private final String cdr3
    private final int sample

    IlmnInfo(int tile, int x, int y, int qSum, String cdr3, int sample) {
        this.tile = tile
        this.x = x
        this.y = y
        this.qSum = qSum
        this.cdr3 = cdr3
        this.sample = sample
    }

    boolean sameTile(IlmnInfo other) {
        this.tile == other.tile
    }

    boolean sameCdr3(IlmnInfo other) {
        this.cdr3 == other.cdr3
    }

    boolean sameSample(IlmnInfo other) {
        this.sample == other.sample
    }

    double dist(IlmnInfo other) {
        int dx = this.x - other.x, dy = this.y - other.y
        Math.sqrt(dx * dx + dy * dy)
    }
}

def umiToInfo = new HashMap<String, List<IlmnInfo>>()

// read in data

int nReads = 0
Double umiSz = null

inputFileNames.eachWithIndex { fileName, ind ->
    println "[${new Date()}] Reading in data from $fileName"

    def reader = new BufferedReader(new InputStreamReader(fileName.endsWith(".gz") ?
            new GZIPInputStream(new FileInputStream(fileName)) :
            new FileInputStream(fileName)))

    def header
    while ((header = reader.readLine()) != null) {
        reader.readLine()
        reader.readLine()
        reader.readLine()

        def splitHeader = header.split(" ")

        // parse coordinates & tile from header

        def coordEntry = splitHeader[0]
        def (tile, x, y) = coordEntry.split(":")[-3..-1].collect { it.toInteger() }

        // parse UMI data from header
        if (splitHeader.length < 3 || !splitHeader[2].startsWith("UMI")) {
            println "Apparently no UMI header in sample.. skipping file"
            break
        }

        def umiEntry = splitHeader[2]
        String umi = umiEntry.split(":")[1] // quality can contain :
        int qSum = 0
        for (int i = umiEntry.length() - umi.length(); i < umiEntry.length(); i++)
            qSum += (int) umiEntry.charAt(i)

        qSum -= 33 * umi.length()

        if (!umiSz)
            umiSz = (double) umi.length()

        // parse CDR3 data from header
        if (splitHeader.length < 4 || !splitHeader[3].startsWith("CDR3")) {
            println "Apparently no CDR3 header in sample.. skipping file"
            break
        }

        def cdr3Entry = splitHeader[3]
        String cdr3 = cdr3Entry.split(":")[1]

        def infoList = umiToInfo[umi]
        if (infoList == null)
            umiToInfo.put(umi, infoList = new ArrayList<IlmnInfo>())

        infoList.add(new IlmnInfo(tile, x, y, qSum, cdr3, ind))

        if (++nReads % 500000 == 0)
            println "[${new Date()}] ${nReads} processed, ${umiToInfo.size()} total UMIs"
    }
}

// Copy to array for random access

println "[${new Date()}] Flattening the map to array for random access"

int n = umiToInfo.size(), ind = 0

def infoArr = new List<IlmnInfo>[n]

umiToInfo.values().each {
    if (it.size() > 1)       // non-trivial MIGs
        infoArr[ind++] = it
    else
        n--
}

// Main cycle. Count stats for umis

println "[${new Date()}] Started parallel processing"

int nPairsToTest = 10_000_000, nRounds = 20, nPairsCapacity = 1_000_000
def THREADS = Runtime.getRuntime().availableProcessors()

def rnd = new Random()

def tileMatch = new AtomicIntegerArray[2]
processed = new AtomicInteger()

tileMatch[0] = new AtomicIntegerArray(4)
tileMatch[1] = new AtomicIntegerArray(4)

def sets = [0: "diff_sample\tdiff_cdr3",
            1: "same_sample\tdiff_cdr3",
            2: "diff_sample\tsame_cdr3",
            3: "same_sample\tsame_cdr3"]

def getSetId = { IlmnInfo info1, IlmnInfo info2 ->
    boolean s = info1.sameSample(info2), c = info1.sameCdr3(info2)
    s ? (c ? 3 : 1) : (c ? 2 : 0)
}

// ind1 - same/diff umi | ind2 - set based on sample/ cdr3 | pair index | dist/qual
def stats = new double[2][4][nPairsCapacity][2]

(0..<nRounds).each {
    println "Round $it"
    GParsPool.withPool THREADS, {
        (0..<nPairsToTest).eachParallel {
            int r1 = rnd.nextInt(n), r2 = rnd.nextInt(n), r3

            while (r1 == r2)
                r2 = rnd.nextInt(n)

            def randomUmiInfoList1 = infoArr[r1],
                randomUmiInfoList2 = infoArr[r2]

            r1 = rnd.nextInt(randomUmiInfoList1.size())
            r2 = rnd.nextInt(randomUmiInfoList1.size())
            r3 = rnd.nextInt(randomUmiInfoList2.size())

            while (r1 == r2)
                r2 = rnd.nextInt(randomUmiInfoList1.size())

            def info1 = randomUmiInfoList1[r1],
                info2 = randomUmiInfoList1[r2], info3 = randomUmiInfoList2[r3]

            // same umi
            if (info1.sameTile(info2)) {
                def setId = getSetId(info1, info2)
                def pairIndex = tileMatch[0].incrementAndGet(setId)

                if (pairIndex < nPairsCapacity) {
                    stats[0][setId][pairIndex][0] = info1.dist(info2)
                    stats[0][setId][pairIndex][1] = Math.min(info1.qSum, info3.qSum) / umiSz
                }
            }

            // diff umi
            if (info1.sameTile(info3)) {
                def setId = getSetId(info1, info3)
                def pairIndex = tileMatch[1].incrementAndGet(setId)

                if (pairIndex < nPairsCapacity) {
                    stats[1][setId][pairIndex][0] = info1.dist(info3)
                    stats[1][setId][pairIndex][1] = Math.min(info1.qSum, info3.qSum) / umiSz
                    tileMatch[1].incrementAndGet(setId)
                }
            }

            int p
            if ((p = processed.incrementAndGet()) % 1000000 == 0)
                println "[${new Date()}] Processed $p of ${nRounds * nPairsToTest} read pairs."
        }
    }
}

// Writing stats

println "[${new Date()}] Writing output"

new File("${outputPrefix}_ilmn_dupl_summary.txt").withPrintWriter { pw ->
    pw.println("umi\t" + sets.values().collect { it.replace('\t', "|") }.join("\t"))
    pw.println("same\t" + (0..3).collect() { tileMatch[0].get(it) }.join("\t"))
    pw.println("diff\t" + (0..3).collect() { tileMatch[1].get(it) }.join("\t"))
    pw.println("total_pairs\t${nRounds * nPairsToTest}")
}

new File("${outputPrefix}_ilmn_dupl_stats.txt").withPrintWriter { pw ->
    pw.println("dist\tqual\tumi\tsample\tcdr3")

    [0: "same_umi", 1: "diff_umi"].each { umi ->
        sets.each { set ->
            stats[umi.key][set.key].each { double[] dq ->
                if (dq[1] > 0 && dq[0] > 0) { // the latter is sampling collision
                    pw.println(dq[0] + "\t" + dq[1] + "\t" + umi.value + "\t" + set.value)
                }
            }
        }
    }
}

println "[${new Date()}] Done"