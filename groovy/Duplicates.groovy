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

package com.milaboratory.migec

@Grab(group = 'org.codehaus.gpars', module = 'gpars', version = '1.2.1')

import groovyx.gpars.GParsPool

import java.util.concurrent.atomic.AtomicInteger
import java.util.zip.GZIPInputStream

/**
 * Input arguments:
 *
 * .FASTQ[.gz] file name list
 * Those should come from migec/Checkout de-multiplexing utility, having UMI:NNN:QQQ in their header
 * where NNN is the UMI sequence and QQQ is its quality Phred score string
 *
 * Output file prefix (coming last)
 */

def inputFileNames = args[0..-2], outputPrefix = args[-1], sameBarcodeSpaceFlag = false

// Core data structure

class IlmnInfo {
    final int tile, x, y, qSum

    IlmnInfo(int tile, int x, int y, int qSum) {
        this.tile = tile
        this.x = x
        this.y = y
        this.qSum = qSum
    }

    boolean sameTile(IlmnInfo other) {
        this.tile == other.tile
    }

    double dist(IlmnInfo other) {
        int dx = this.x - other.x, dy = this.y - other.y
        Math.sqrt(dx * dx + dy * dy)
    }
}

def umiToInfo = new HashMap<String, List<IlmnInfo>>()

// read in data

inputFileNames.eachWithIndex { fileName, ind ->
    println "[${new Date()}] Reading in data from $fileName"

    def reader = new BufferedReader(new InputStreamReader(fileName.endsWith(".gz") ?
            new GZIPInputStream(new FileInputStream(fileName)) :
            new FileInputStream(fileName)))
    int nReads = 0
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

        def umiEntry = splitHeader[2]
        String umi = umiEntry.split(":")[1] // quality can contain :
        if (!sameBarcodeSpaceFlag)
            umi += ind // to distinguish different samples
        int qSum = 0
        for (int i = umiEntry.length() - umi.length(); i < umiEntry.length(); i++)
            qSum += (int) umiEntry.charAt(i)

        qSum -= 33 * umi.length()

        def infoList = umiToInfo[umi]
        if (infoList == null)
            umiToInfo.put(umi, infoList = new ArrayList<IlmnInfo>())

        infoList.add(new IlmnInfo(tile, x, y, qSum))

        if (++nReads % 500000 == 0)
            println "[${new Date()}] ${nReads} processed, ${umiToInfo.size()} total UMIs"
    }
}

def umiSz = (double) umiToInfo.keySet().iterator().next().length()

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

int nPairsToTest = 10000000
def THREADS = Runtime.getRuntime().availableProcessors()

def rnd = new Random()

def sameTileCount = new AtomicInteger(), bgSameTileCount = new AtomicInteger(), processed = new AtomicInteger()
def stats = new double[nPairsToTest][4]

GParsPool.withPool THREADS, {
    (0..<nPairsToTest).eachParallel { int pairIndex ->

        def randomUmiInfoList1 = infoArr[rnd.nextInt(n)],
            randomUmiInfoList2 = infoArr[rnd.nextInt(n)]

        int r1 = rnd.nextInt(randomUmiInfoList1.size()),
            r2 = rnd.nextInt(randomUmiInfoList1.size()),
            r3 = rnd.nextInt(randomUmiInfoList2.size())

        while (r1 == r2)
            r2 = rnd.nextInt(randomUmiInfoList1.size())

        def info1 = randomUmiInfoList1[r1],
            info2 = randomUmiInfoList1[r2], info3 = randomUmiInfoList2[r3]

        boolean sameTile = info1.sameTile(info2),
                bgSameTile = info1.sameTile(info3)
        stats[pairIndex][0] = sameTile ? info1.dist(info2) : Double.NaN
        stats[pairIndex][1] = bgSameTile ? info1.dist(info3) : Double.NaN
        stats[pairIndex][2] = sameTile ? Math.min(info1.qSum, info2.qSum) / umiSz : Double.NaN // by umi length
        stats[pairIndex][3] = bgSameTile ? Math.min(info1.qSum, info3.qSum) / umiSz : Double.NaN // by umi length

        if (sameTile)
            sameTileCount.incrementAndGet()

        if (bgSameTile)
            bgSameTileCount.incrementAndGet()

        int p
        if ((p = processed.incrementAndGet()) % 500000 == 0)
            println "[${new Date()}] Processed $p of $nPairsToTest read pairs. " +
                    "Tile matched (approx) ${sameTileCount.get()} for same UMI and ${bgSameTileCount.get()} for random."
    }
}

// Writing stats

println "[${new Date()}] Writing output"

new File("${outputPrefix}_ilmn_dupl_summary.txt").withPrintWriter { pw ->
    pw.println("same_tile\t${sameTileCount.get()}")
    pw.println("bg_same_tile\t${bgSameTileCount.get()}")
    pw.println("total_pairs\t${nPairsToTest}")
}

new File("${outputPrefix}_ilmn_dupl_stats.txt").withPrintWriter { pw ->
    pw.println("dist\tbg.dist\tqual\tbg.qual")

    for (int i = 0; i < nPairsToTest; i++) {
        double[] arr = stats[i]
        if (!arr[0].isNaN() || !arr[1].isNaN())
            pw.println(arr.collect { it.isNaN() ? "NA" : it }.join("\t"))
    }
}

println "[${new Date()}] Done"