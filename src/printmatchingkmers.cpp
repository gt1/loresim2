/*
    loresim2
    Copyright (C) 2015-2016 German Tischler

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
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/util/PrefixSums.hpp>

struct AlignmentReaderWrapper
{
	libmaus2::bambam::BamDecoder & dec_a;
	int64_t id;

	AlignmentReaderWrapper(libmaus2::bambam::BamDecoder & rdec_a)
	: dec_a(rdec_a), id(-1)
	{

	}

	bool readAlignment()
	{
		bool const ok = dec_a.readAlignment();
		if ( ok )
			++id;
		return ok;
	}
};

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		std::string const fn_a = arg[0];
		std::set<std::string> S;
		uint64_t found = 0;
		for ( uint64_t i = 1 ; i < arg.size(); ++i )
			S.insert(arg[i]);
		uint64_t const k = arg.uniqueArgPresent("k") ? arg.getUnsignedNumericArg<uint64_t>("k") : 20;

		libmaus2::bambam::BamDecoder dec_a(fn_a);
		AlignmentReaderWrapper ARW(dec_a);
		libmaus2::bambam::BamAlignment const & algn = dec_a.getAlignment();
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Acigop;
		libmaus2::lcs::AlignmentTraceContainer ATC;

		while ( ((!S.size()) || (found < S.size())) && ARW.readAlignment() )
		{
			if ( algn.isMapped() )
			{
				bool printedname = false;

				if ( (!S.size()) || (S.find(algn.getName()) != S.end()) )
				{
					uint32_t const numcig = algn.getCigarOperations(Acigop);
					libmaus2::bambam::CigarStringParser::cigarToTrace(Acigop.begin(),Acigop.begin()+numcig, ATC, true);
					uint64_t const fclip = libmaus2::bambam::BamAlignmentDecoderBase::getFrontClipping(algn.D.begin());
					uint64_t const seqlen = dec_a.getHeader().getRefIDLength(algn.getRefID());
					uint64_t const readlength = algn.getLseq();

					uint64_t apos = algn.getPos();
					uint64_t bpos = fclip;
					libmaus2::lcs::AlignmentTraceContainer::step_type const * tc = ATC.ta;
					while ( tc != ATC.te )
					{
						switch ( *tc )
						{
							case libmaus2::lcs::AlignmentTraceContainer::STEP_INS:
								bpos++;
								tc++;
								break;
							case libmaus2::lcs::AlignmentTraceContainer::STEP_DEL:
								apos++;
								tc++;
								break;
							case libmaus2::lcs::AlignmentTraceContainer::STEP_MISMATCH:
								apos++;
								bpos++;
								tc++;
								break;
							case libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH:
							{
								uint64_t cnt = 0;

								while ( tc != ATC.te && *tc == libmaus2::lcs::AlignmentTraceContainer::STEP_MATCH )
								{
									cnt += 1;
									apos += 1;
									bpos += 1;
									tc++;
								}

								if ( cnt >= k )
								{
									if ( ! printedname )
									{
										std::cout
											<< "refid=" << algn.getRefID()
											<< "\trefname=" << dec_a.getHeader().getRefIDName(algn.getRefID())
											<< "\trefseqlen=" << seqlen
											<< "\trefalnlen=" << algn.getReferenceLength()
											<< "\treadid=" << ARW.id
											<< "\treadname=" << algn.getName()
											<< "\treadlen=" << readlength
											<< "\tdir=" << (algn.isReverse()?"reco":"forw")
											<< std::endl;
										if ( algn.getAuxString("er") )
											std::cout << algn.getAuxString("er") << std::endl;
									}
									std::cout << cnt << "\t" << apos-cnt << "\t" << bpos-cnt << std::endl;
									printedname = true;
								}
								break;
							}
							default:
								break;
						}
					}

					#if 0
					std::vector < std::pair<uint64_t,uint64_t> > const off =
						ATC.getKMatchOffsets(k, algn.getPos(), fclip);

                                        int64_t prev = -1;
					for ( uint64_t i = 0; i < off.size(); ++i )
					{
						std::cout << algn.getName() << "\t" << fclip << "\t" << algn.getRefID() << "\t" << off[i].first << "\t" << off[i].second << "\t" << off[i].second-prev << std::endl;
						prev = off[i].second;
					}
					#endif

					found += 1;
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
