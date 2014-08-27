/* Copyright (c) 2010-2012 Sebastian Bauer
 *
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted (subject to the limitations in the
 * disclaimer below) provided that the following conditions are met:
 *
 * * Redistributions of source code must retain the above copyright
 *   notice, this list of conditions and the following disclaimer.
 *
 * * Redistributions in binary form must reproduce the above copyright
 *   notice, this list of conditions and the following disclaimer in the
 *   documentation and/or other materials provided with the
 *   distribution.
 *
 * * Neither the name of Sebastian Bauer nor the names of its
 *   contributors may be used to endorse or promote products derived
 *   from this software without specific prior written permission.
 *
 * NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
 * GRANTED BY THIS LICENSE.  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT
 * HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
 * WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
 * BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
 * WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE
 * OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN
 * IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package sonumina.boqa;

import java.io.IOException;

import ontologizer.GlobalPreferences;
import ontologizer.OntologizerThreadGroups;
import ontologizer.benchmark.Datafiles;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.GnuParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import sonumina.boqa.benchmark.Benchmark;
import sonumina.boqa.calculation.BOQA;

/**
 * Main entry point of the BOQA benchmark.
 * 
 * @author Sebastian Bauer
 */
public class BOQABenchmark
{
	static private String ontologyPath;
	static private String annotationPath;

	static private double ALPHA = 0.002;
	static private double BETA = 0.1;
	static private int MAX_TERMS = -1;
	static private int SAMPLES_PER_ITEM = 5;
	static private boolean CONSIDER_FREQUENCIES_ONLY = false;
	static private int SIZE_OF_SCORE_DISTRIBUTION = 250000;
	static private String RESULT_BASE_NAME = "benchmark";
	
	/**
	 * Parses the command line and returns a corresponding
	 * BOQA object. 
	 * 
	 * @param args
	 */
	public static BOQA parseCommandLine(String [] args)
	{
	   Options opt = new Options();
	   opt.addOption("o", "ontology", true, "Path or URL to the ontology file.");
	   opt.addOption("a", "annotations", true, "Path or URL to files containing annotations.");
	   opt.addOption("c", "considerFreqOnly", false, "If specified, only items with frequencies are considered.");
	   opt.addOption("m", "maxTerms", true, "Defines the maximal number of terms a random query can have. Default is " + MAX_TERMS);
	   opt.addOption("s", "samplesPerItem", true, "Define the number of samples per item. Defaults to " + SAMPLES_PER_ITEM + ".");
	   opt.addOption("r", "resultBaseName", true, "Defines the base name of the result files that are created during the benchmark. Defaults to \"" + RESULT_BASE_NAME + "\".");
	   opt.addOption(null, "alpha", true, "Specifies alpha (false-positive rate) during simulation. Default is " + ALPHA + ".");
	   opt.addOption(null, "beta", true, "Specifies beta (false-negative rate) during simulation. Default is " + BETA + ".");
	   opt.addOption(null, "sizeOfScoreDistribution", true, "Specifies the size of the score distribution. Default is " + SIZE_OF_SCORE_DISTRIBUTION + ".");
	   opt.addOption("h", "help", false, "Shows this help");

	   BOQA boqa = new BOQA();

	   try
	   {
		   GnuParser parser = new GnuParser();
		   CommandLine cl;
		   cl = parser.parse(opt, args);

		   if (cl.hasOption('h'))
		   {
			   HelpFormatter f = new HelpFormatter();
			   f.printHelp(BOQABenchmark.class.getName(), opt);
			   System.exit(0);
		   }
		   
		   if (cl.hasOption('m'))
			   MAX_TERMS = Integer.parseInt(cl.getOptionValue('m'));

		   if (cl.hasOption('c'))
			   CONSIDER_FREQUENCIES_ONLY = true;

		   if (cl.hasOption('s'))
			   SAMPLES_PER_ITEM = Integer.parseInt(cl.getOptionValue('s'));

		   SIZE_OF_SCORE_DISTRIBUTION = Integer.parseInt(cl.getOptionValue("sizeOfScoreDistribution", "250000"));
		   RESULT_BASE_NAME = cl.getOptionValue('r', RESULT_BASE_NAME);
		   
		   if (cl.hasOption("alpha"))
			   ALPHA = Double.parseDouble(cl.getOptionValue("alpha"));

		   if (cl.hasOption("beta"))
			   BETA = Double.parseDouble(cl.getOptionValue("beta"));

		   ontologyPath = cl.getOptionValue('o',ontologyPath);
		   annotationPath = cl.getOptionValue('a', annotationPath);
		   
		   boqa.setSimulationAlpha(ALPHA);
		   boqa.setSimulationBeta(BETA);
		   boqa.setConsiderFrequenciesOnly(CONSIDER_FREQUENCIES_ONLY);
		   boqa.setSimulationMaxTerms(MAX_TERMS);
		   if (MAX_TERMS != -1)
			   boqa.setMaxQuerySizeForCachedDistribution(MAX_TERMS);
		   boqa.setSizeOfScoreDistribution(SIZE_OF_SCORE_DISTRIBUTION);
	   } catch (ParseException e)
	   {
		   System.err.println("Faield to parse commandline: " + e.getLocalizedMessage());
		   System.exit(1);
	   }
	   return boqa;
	}

	/**
	 * The main entry.
	 *  
	 * @param args
	 * @throws IOException 
	 * @throws InterruptedException 
	 */
	public static void main(String[] args) throws InterruptedException, IOException
	{
		BOQA boqa = parseCommandLine(args);

		boqa.setPrecalculateJaccard(true);

		GlobalPreferences.setProxyPort(888);
		GlobalPreferences.setProxyHost("realproxy.charite.de");

		Datafiles df = new Datafiles(ontologyPath,annotationPath);
		boqa.setup(df.graph, df.assoc);
		
		Benchmark benchmark = new Benchmark();
		benchmark.setSamplesPerItem(SAMPLES_PER_ITEM);
		benchmark.setResultBaseName(RESULT_BASE_NAME);
		benchmark.benchmark(boqa);

		OntologizerThreadGroups.workerThreadGroup.interrupt();
	}
}
