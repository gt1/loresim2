/*
    loresim
    Copyright (C) 2015 German Tischler

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <libmaus2/bambam/BamHeader.hpp>
#include <libmaus2/bambam/BamBlockWriterBaseFactory.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/random/DNABaseNoiseSpiker.hpp>
#include <libmaus2/util/ArgInfo.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/fastx/FastAIndexGenerator.hpp>
#include "Intv.hpp"

std::string getRandom(uint64_t const randlen)
{
	std::ostringstream o;

	for ( uint64_t i = 0; i < randlen; ++i )
	{
        	switch ( libmaus2::random::Random::rand8() % 4 )
                {
                	case 0: o.put('A'); break;
                	case 1: o.put('C'); break;
                	case 2: o.put('G'); break;
                	case 3: o.put('T'); break;
		}
	}

	return o.str();
}

int main(int argc, char ** argv)
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);
		uint64_t randseed = arginfo.getValueUnsignedNumeric<uint64_t>("randomseed",time(0));
		libmaus2::random::Random::setup(randseed);

		libmaus2::bambam::BamSeqEncodeTable const seqenc;
		::libmaus2::fastx::UCharBuffer UB;

		// substitution error fraction
		double substrate = arginfo.getValue<double>("substrate",.2);
		// deletion error fraction
		double delrate = arginfo.getValue<double>("delrate",.3);
		// insertion error fraction
		double insrate = arginfo.getValue<double>("insrate",.5);
		// insertion homopolymer error fraction
		double inshomopolrate = arginfo.getValue<double>("inshomopolrate",.0);

		// low error rate
		double eratelow        = arginfo.getValue<double>("eratelow",.15);
		// high error rate
		double eratehigh       = arginfo.getValue<double>("eratehigh",.25);
		// std dev for low error rate
		double eratelowstddev  = arginfo.getValue<double>("eratelowstddev",0.03);
		// std dev for high error rate
		double eratehighstddev = arginfo.getValue<double>("eratehighstddev",0.04);

		// drop rate
		double const droprate = arginfo.getValue<double>("droprate",0.00);
		// number of traversals per input sequence
		uint64_t const numtraversals = arginfo.getValueUnsignedNumeric<uint64_t>("numtraversals",1);
		// average read length
		double const readlenavg = arginfo.getValue<double>("readlenavg", 15000);
		// read length std dev
		double const readlenstddev = arginfo.getValue<double>("readlenstddev", 3000);
		// probability to stay in low error rate mode
		double const keeplowstate = arginfo.getValue<double>("keeplowstate", 0.9998);
		// probability to stay in high error rate mode
		double const keephighstate = arginfo.getValue<double>("keeplowstate", 0.995);
		// probability to start in low error rate mode
		double const startlowprob = arginfo.getValue<double>("startlowprob", 0.7);
		// number of random bases appended at front and back
		uint64_t const randlen = arginfo.getValueUnsignedNumeric<uint64_t>("randlen",500);
		bool const placerandom = arginfo.getValueUnsignedNumeric<uint64_t>("placerandom",0);
		bool const nthres = arginfo.getValueUnsignedNumeric<uint64_t>("nthres",0);

		uint64_t const minlen = arginfo.getValueUnsignedNumeric<uint64_t>("minlen",0);

		// noise spiker object
		libmaus2::random::DNABaseNoiseSpiker DBNS(substrate, delrate, insrate, inshomopolrate, eratelow, eratehigh, eratelowstddev, eratehighstddev, keeplowstate, keephighstate, startlowprob);

		// input reference name (FastA format)
		std::string const reffn = arginfo.restargs.at(0);
		// output fasta file name
		std::string const fafn = arginfo.restargs.at(1);

		// construct name of FAI file
		std::string const fainame = reffn + ".fai";
		libmaus2::fastx::FastAIndexGenerator::generate(reffn,fainame,1);
		// load FAI index for reference FastA
		libmaus2::fastx::FastAIndex::unique_ptr_type PFAI(libmaus2::fastx::FastAIndex::load(fainame));
		// open output fasta file
		libmaus2::aio::OutputStreamInstance faOSI(fafn);

		// create SAM header
		std::ostringstream samheaderstr;
		samheaderstr << "@HD\tVN:1.5\tSO:unknown\n";
		for ( uint64_t i = 0; i < PFAI->size(); ++i )
		{
			libmaus2::fastx::FastAIndexEntry const & entry = (*PFAI)[i];
			samheaderstr << "@SQ\tSN:" << entry.name << "\tLN:" << entry.length << "\n";
		}

		libmaus2::bambam::BamHeader header(samheaderstr.str());

		// bam writer
		libmaus2::bambam::BamBlockWriterBase::unique_ptr_type writer(
			libmaus2::bambam::BamBlockWriterBaseFactory::construct(header, arginfo)
		);

		// id of run
		uint64_t const runid = 0;
		// current well id
		uint64_t wellid = 0;

		libmaus2::aio::InputStreamInstance ISI(reffn);
		libmaus2::fastx::StreamFastAReaderWrapper SFAR(ISI);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;

		for ( uint64_t seqid = 0; SFAR.getNextPatternUnlocked(pattern); ++seqid )
		{
			std::cerr << "[V] " << pattern.sid << std::endl;

			std::ostringstream seqstr;
			seqstr << pattern.spattern;

			std::string const seq = seqstr.str();

			for ( uint64_t trav = 0; trav < numtraversals; ++trav )
			{
				std::cerr << "\t[V] traversal " << trav << std::endl;

				uint64_t pp = 0;
				typedef std::pair<uint64_t,uint64_t> upair;
				std::vector<upair> poslenvec;

				if ( placerandom )
					while ( pp < seq.size() )
					{
						int64_t len = libmaus2::random::GaussianRandom::random(readlenstddev,readlenavg);
						len = std::min(len,static_cast<int64_t>(seq.size()));
						uint64_t p = libmaus2::random::Random::rand64() % (seq.size()-len+1);
						poslenvec.push_back(upair(p,len));
						pp += len;
					}
				else
					while ( pp < seq.size() )
					{
						int64_t len = libmaus2::random::GaussianRandom::random(readlenstddev,readlenavg);

						len = std::min(len,static_cast<int64_t>(seq.size()-pp));

						poslenvec.push_back(upair(pp,len));

						pp += len;
					}

				for ( uint64_t z = 0; z < poslenvec.size(); ++z )
				{
					if ( libmaus2::random::UniformUnitRandom::uniformUnitRandom() < droprate )
						continue;
					uint64_t const len = poslenvec[z].second;
					if ( len < minlen )
						continue;

					uint64_t const pos = poslenvec[z].first;
					bool const strand = libmaus2::random::Random::rand8() & 1;
					bool rc = !strand;
					std::string sub = seq.substr(pos,len);

					uint64_t ncnt = 0;
					for ( uint64_t i = 0; i < sub.size(); ++i )
						if ( sub[i] == 'N' )
							ncnt++;
					if ( ncnt > nthres )
						continue;

					libmaus2::random::DNABaseNoiseSpiker::ErrorStats E;
					std::pair<std::string,std::string> const P = DBNS.modifyAndComment(sub,&E);

					// std::cerr << P.second << std::endl;

					std::string errintv = P.second;
					if ( errintv.find("ERRINTV=[[") != std::string::npos )
					{
						errintv = errintv.substr(errintv.find("ERRINTV=[[" )+ strlen("ERRINTV=[["));
						if ( errintv.find("]]") != std::string::npos )
						{
							errintv = errintv.substr(0,errintv.find("]]"));
							#if 0
							std::deque<std::string> Vintv(libmaus2::util::stringFunctions::tokenize<std::string>(errintv,";"));
							for ( uint64_t i = 0; i < Vintv.size(); ++i )
								std::cerr << Vintv[i] << std::endl;
							#endif
						}
						else
						{
							errintv = std::string();
						}
					}
					else
					{
						errintv = std::string();
					}

					std::string mod = P.first;

					std::string const r = getRandom(randlen) + mod + getRandom(randlen);

					std::string cig = P.second;
					if ( cig.find("CIGAR=[") != std::string::npos )
						cig = cig.substr(cig.find("CIGAR=[")+strlen("CIGAR=["));
					if ( cig.find("]") != std::string::npos )
						cig = cig.substr(0,cig.find("]"));

					std::vector<libmaus2::bambam::cigar_operation> Vcigop = libmaus2::bambam::CigarStringParser::parseCigarString(cig);

					uint64_t delshift = 0;
					for (
						uint64_t ind = 0;
						ind < Vcigop.size() &&
						Vcigop[ind].first != libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH &&
						Vcigop[ind].first != libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL &&
						Vcigop[ind].first != libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF
						; ++ind
					)
					{
						if ( Vcigop[ind].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL )
							delshift += Vcigop[ind].second;
					}

					std::ostringstream padcigstr;
					if ( randlen )
						padcigstr << randlen << "S";
					padcigstr << cig;
					if ( randlen )
						padcigstr << randlen << "S";
					cig = padcigstr.str();

					std::ostringstream namestr;
					std::string const rfasta = strand ? r : libmaus2::fastx::reverseComplementUnmapped(r);
					namestr << 'L' << runid << '/' << (wellid++) << '/' << 0 << '_' << rfasta.size();
					std::string const name = namestr.str();

					UB.reset();
					std::vector<char> Vqual(r.size(),'H');

					try
					{
						libmaus2::bambam::BamAlignmentEncoderBase::encodeAlignment
						(
							UB,
							seqenc,
							name,
							seqid,
							pos + delshift,
							255,
							rc ? libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_FREVERSE : 0,
							cig,
							-1,
							-1,
							0,
							r,
							std::string(Vqual.begin(),Vqual.end()),
							33,true
						);
					}
					catch(std::exception const & lme)
					{
						std::cerr << lme.what() << std::endl;
						continue;
					}

					::libmaus2::bambam::MdStringComputationContext mdnmcontext;
					libmaus2::bambam::BamAlignmentDecoderBase::calculateMd(UB.buffer,UB.length,mdnmcontext,seq.begin()+pos /* itref */);

					libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(UB,"MD",mdnmcontext.md.get());
					libmaus2::bambam::BamAlignmentEncoderBase::putAuxNumber(UB,"NM",'i',mdnmcontext.nm);

					if ( errintv.size() )
					{
						libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(
							UB,
							"er",
							errintv.c_str()
						);

						std::vector<Intv> Vintv = Intv::parse(errintv);
						std::vector<Intv> Rintv = Intv::computeRIntv(Vintv,UB.buffer);

						std::ostringstream ostr;
						for ( uint64_t i = 0; i < Rintv.size(); ++i )
							ostr << Rintv[i] << ";";

						libmaus2::bambam::BamAlignmentEncoderBase::putAuxString(
							UB,
							"ee",
							ostr.str().c_str()
						);
					}

					faOSI << '>' << name << " RQ=0.851" << " POS=[" << pos << "] REFID=[" << seqid << "] RC=[" << rc << "]\n";
					for ( uint64_t i = 0; i < rfasta.size(); )
					{
						uint64_t toprint = std::min(static_cast<int>(rfasta.size()-i),static_cast<int>(80));
						faOSI.write(rfasta.c_str()+i,toprint);
						i += toprint;
						faOSI.put('\n');
					}
					writer->writeBamBlock(UB.buffer,UB.length);
				}
			}
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
