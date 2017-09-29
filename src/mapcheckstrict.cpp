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
#include <libmaus2/lcs/NNPLocalAligner.hpp>
#include <libmaus2/util/ArgParser.hpp>
#include <libmaus2/bambam/BamDecoder.hpp>
#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/dazzler/align/OverlapIndexer.hpp>
#include <libmaus2/fastx/FastAReader.hpp>
#include <libmaus2/util/PrefixSums.hpp>

#include "Intv.hpp"

struct OverlapBComparator
{
	bool operator()(libmaus2::dazzler::align::Overlap const & A, libmaus2::dazzler::align::Overlap const & B) const
	{
		return A.bread < B.bread;
	}
};

std::string getPaddedRead(libmaus2::fastx::FastAReader::pattern_type const & pattern, bool const rc)
{
	std::ostringstream ostr;

	ostr.put(4);
	for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
		ostr.put(libmaus2::fastx::mapChar(pattern.spattern[i]));
	ostr.put(4);

	if ( rc )
	{
		ostr.put(4);
		for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
		{
			uint64_t const rpos = pattern.spattern.size()-i-1;
			char const base = pattern.spattern[rpos];
			char const premappedbase = libmaus2::fastx::mapChar(base);
			char const mappedbase = premappedbase < 4 ? premappedbase : 0;
			char const invbase = libmaus2::fastx::invertN(mappedbase);
			ostr.put(invbase);
		}
		ostr.put(4);
	}

	return ostr.str();
}

void loadfasta(std::string const rnfa, std::map<std::string,uint64_t> & readnametoid, std::vector<uint64_t> & readoff, std::string & readdata)
{
	std::cerr << "[V] loading " << rnfa << "...";
	libmaus2::fastx::FastAReader FR(rnfa);
	libmaus2::fastx::FastAReader::pattern_type pattern;
	uint64_t id = 0;
	std::ostringstream readostr;
	while ( FR.getNextPatternUnlocked(pattern) )
	{
		std::string const shortid = pattern.getShortStringId();
		readnametoid[shortid] = id++;
		std::string const padded = getPaddedRead(pattern,true);
		readostr << padded;
		readoff.push_back(padded.size());
	}

	uint64_t const s = libmaus2::util::PrefixSums::prefixSums(readoff.begin(),readoff.end());
	readoff.push_back(s);
	readdata = readostr.str();
	std::cerr << "done." << std::endl;
}

#include <libmaus2/lcs/NP.hpp>

std::ostream & operator<<(std::ostream& out, std::pair<uint64_t,uint64_t> const & P)
{
	out << "P(" << P.first << "," << P.second << ")";
	return out;
}

struct BamPeeker
{
	libmaus2::bambam::BamDecoder & bamdec;
	libmaus2::bambam::BamAlignment const & algn;
	bool have;

	void consume()
	{
		have = bamdec.readAlignment();
	}

	BamPeeker(libmaus2::bambam::BamDecoder & rbamdec) : bamdec(rbamdec), algn(bamdec.getAlignment())
	{
		consume();
	}

	bool peek()
	{
		return have;
	}

	std::string getName() const
	{
		return algn.getName();
	}

	bool getNextRange(std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> & V)
	{
		V.resize(0);

		if ( peek() )
		{
			std::string const refname = getName();

			while ( peek() && getName() == refname )
			{
				V.push_back(algn.sclone());
				consume();
			}

			return true;
		}
		else
		{
			return false;
		}
	}
};

double handleBlock(std::vector<libmaus2::dazzler::align::Overlap> & VOVL, libmaus2::bambam::BamDecoder & dec_a, int64_t & rid, int64_t & rlen)
{
	int64_t const bread = VOVL.front().bread;
	while ( rid < bread )
	{
		rid += 1;
		bool const ok = dec_a.readAlignment();
		assert ( ok );
		rlen = dec_a.getAlignment().getLseq();
		rlen -= (
			dec_a.getAlignment().getFrontSoftClipping()+
			dec_a.getAlignment().getBackSoftClipping()
		);
	}
	assert ( rid == bread );
	std::vector< libmaus2::math::IntegerInterval<int64_t> > Vintv;
	for ( uint64_t i = 0; i < VOVL.size(); ++i )
	{
		libmaus2::math::IntegerInterval<int64_t> intv(VOVL[i].path.bbpos,VOVL[i].path.bepos-1);
		Vintv.push_back(intv);
	}
	Vintv = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vintv);
	int64_t diam = 0;
	for ( uint64_t i = 0; i < Vintv.size(); ++i )
		diam += Vintv[i].diameter();

	double const frac = static_cast<double>(diam) / rlen;
	//std::cerr << frac << std::endl;

	return frac;
}

struct FastAGetter
{
	libmaus2::fastx::FastAReader::unique_ptr_type p;
	libmaus2::fastx::FastAReader::pattern_type pattern;

	FastAGetter(std::string const & fn) : p(fn.size() ? new libmaus2::fastx::FastAReader(fn) : NULL)
	{

	}

	void getNext()
	{
		if ( p )
			p->getNextPatternUnlocked(pattern);
	}
};

#include <libmaus2/bambam/BamWriter.hpp>

struct MissingOutput
{
	typedef MissingOutput this_type;
	typedef libmaus2::util::unique_ptr<this_type>::type unique_ptr_type;

	std::string const bamfilename;
	libmaus2::bambam::BamWriter::unique_ptr_type BW;
	std::string const fastafilename;
	libmaus2::aio::OutputStreamInstance::unique_ptr_type OSI;

	MissingOutput(std::string const & filename, libmaus2::bambam::BamHeader const & header, bool const open)
	:
		bamfilename(filename+".bam"),
		fastafilename(filename+".fasta")
	{
		if ( open )
		{
			libmaus2::bambam::BamWriter::unique_ptr_type tBW(new libmaus2::bambam::BamWriter(bamfilename,header));
			BW = UNIQUE_PTR_MOVE(tBW);

			libmaus2::aio::OutputStreamInstance::unique_ptr_type tOSI(new libmaus2::aio::OutputStreamInstance(fastafilename));
			OSI = UNIQUE_PTR_MOVE(tOSI);
		}
	}

	void put(libmaus2::bambam::BamAlignment const & algn, libmaus2::fastx::Pattern const & pattern)
	{
		if ( BW )
			BW->writeAlignment(algn);
		if ( OSI )
			(*OSI) << pattern;
	}

	void put(libmaus2::bambam::BamAlignment const & algn)
	{
		if ( BW )
			BW->writeAlignment(algn);
	}
};

#include <libmaus2/util/stringFunctions.hpp>


struct RepeatLine
{
	uint64_t s;
	uint64_t p;
	uint64_t l;
	double e;

	RepeatLine()
	{

	}

	RepeatLine(
		uint64_t rs,
		uint64_t rp,
		uint64_t rl,
		double re
	) : s(rs), p(rp), l(rl), e(re)
	{

	}

	int64_t getFrom() const
	{
		return p;
	}

	int64_t getTo() const
	{
		return p+l;
	}

	libmaus2::math::IntegerInterval<int64_t> getInterval() const
	{
		return libmaus2::math::IntegerInterval<int64_t>(getFrom(),getTo()-1);
	}

	libmaus2::math::IntegerInterval<int64_t> getOverlap(libmaus2::bambam::BamAlignment const & A) const
	{
		return A.getReferenceInterval().intersection(getInterval());
	}

	static uint64_t parseInt(std::string const & s)
	{
		std::istringstream istr(s);
		uint64_t i;
		istr >> i;

		if ( istr.bad() || istr.peek() != std::istream::traits_type::eof() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Cannot parse " << s << " as integer" << std::endl;
			lme.finish();
			throw lme;
		}

		return i;
	}

	static double parseDouble(std::string const & s)
	{
		std::istringstream istr(s);
		double i;
		istr >> i;

		if ( istr.bad() || istr.peek() != std::istream::traits_type::eof() )
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Cannot parse " << s << " as double" << std::endl;
			lme.finish();
			throw lme;
		}

		return i;
	}

	RepeatLine(std::string const & q)
	{
		std::deque<std::string> const T = libmaus2::util::stringFunctions::tokenize<std::string>(q,",");

		if ( T.size() == 4 )
		{
			s = parseInt(T[0]);
			p = parseInt(T[1]);
			l = parseInt(T[2]);
			e = parseDouble(T[3]);
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Cannot parse " << q << " as RepeatLine" << std::endl;
			lme.finish();
			throw lme;
		}
	}
};

std::ostream & operator<<(std::ostream & out, RepeatLine const & RL)
{
	out << "RepeatLine(s=" << RL.s << ",p=" << RL.p << ",l=" << RL.l << ",e=" << RL.e << ")";
	return out;
}


std::vector < RepeatLine > loadRepeatLines(std::string const & fn)
{
	libmaus2::aio::InputStreamInstance ISI(fn);
	std::vector < RepeatLine > VRL;

	while ( ISI )
	{
		std::string line;
		std::getline(ISI,line);
		if ( line.size() )
		{
			RepeatLine RL(line);
			//std::cerr << RL << std::endl;
			VRL.push_back(RL);
		}
	}

	std::cerr << "[V] loaded " << VRL.size() << " repeat lines" << std::endl;

	return VRL;
}

#include <libmaus2/geometry/RangeSet.hpp>
#include <libmaus2/fastx/FastAIndex.hpp>
#include <libmaus2/rank/DNARank.hpp>

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		std::string const fn_a = arg[0];
		std::string const fn_b = arg[1];
		std::string fn_fasta;
		if ( 2 < arg.size() )
			fn_fasta = arg[2];
		FastAGetter FAG(fn_fasta);

		// error fragments we count at all
		double const errthres = arg.uniqueArgPresent("e") ? arg.getParsedArg<double>("e") : 0.15;
		std::string const repmapfn = arg.uniqueArgPresent("R") ? arg["R"] : std::string();
		std::string const indexfn = arg.uniqueArgPresent("I") ? arg["I"] : std::string();
		std::vector <RepeatLine> VRL;
		std::map < uint64_t , libmaus2::geometry::RangeSet<RepeatLine>::shared_ptr_type > remapM;
		uint64_t const mapk = arg.uniqueArgPresent("k") ? arg.getParsedArg<uint64_t>("k") : 20;
		bool const writealgncat = arg.uniqueArgPresent("writealgncat") ? arg.getParsedArg<uint64_t>("writealgncat") : false;

		libmaus2::rank::DNARank::unique_ptr_type Prank;
		if ( indexfn.size() )
		{
			std::cerr << "[V] loading index...";
			libmaus2::rank::DNARank::unique_ptr_type Trank(libmaus2::rank::DNARank::loadFromRunLength(indexfn, 32/*num threads */));
			Prank = UNIQUE_PTR_MOVE(Trank);
			std::cerr << "done\n";
		}

		std::string const fainame = arg.uniqueArgPresent("F") ? arg["F"] : std::string();
                libmaus2::fastx::FastAIndex::unique_ptr_type Pindex;
                if ( fainame.size() )
                {
	                libmaus2::fastx::FastAIndex::unique_ptr_type Tindex(libmaus2::fastx::FastAIndex::load(fainame));
			Pindex = UNIQUE_PTR_MOVE(Tindex);
		}
                // libmaus2::fastx::FastAIndex const & FAI = *Pindex;


		if ( repmapfn.size() )
		{
			VRL = loadRepeatLines(repmapfn);
			uint64_t l = 0;
			while ( l < VRL.size() )
			{
				uint64_t h = l+1;
				while ( h < VRL.size() && VRL[h].s == VRL[l].s )
					++h;

				uint64_t const s = VRL[l].s;

				assert ( Pindex );
				assert ( s < Pindex->size() );

				libmaus2::geometry::RangeSet<RepeatLine>::shared_ptr_type P(
					new libmaus2::geometry::RangeSet<RepeatLine>((*Pindex)[s].length)
				);

				for ( uint64_t i = l; i < h; ++i )
				{
					P->insert(VRL[i]);
				}

				remapM[s] = P;

				l = h;
			}
		}
		// allow this many low error bases to be missed without reporting
		uint64_t const lowerrallow = arg.uniqueArgPresent("L") ? arg.getUnsignedNumericArg<uint64_t>("L") : 150;

		uint64_t verbmod = arg.uniqueArgPresent("V") ? arg.getUnsignedNumericArg<uint64_t>("V") : 1024;
		std::string reference;
		if ( arg.uniqueArgPresent("r") )
			reference = arg["r"];
		std::map<int64_t,std::string> Mref;
		if ( reference.size() )
		{
			libmaus2::fastx::FastAReader F(reference);
			libmaus2::fastx::FastAReader::pattern_type pattern;
			uint64_t id = 0;
			while ( F.getNextPatternUnlocked(pattern) )
			{
				for ( uint64_t i = 0; i < pattern.spattern.size(); ++i )
					pattern.spattern[i] = libmaus2::fastx::mapChar(pattern.spattern[i]);
				Mref[id++] = pattern.spattern;
			}
		}


		#if 0
		std::string const ovl_a = arg[2];
		std::string const reffa = arg[3];
		std::string const rnfa = arg[4];
		#endif

		#if 0
		std::map<std::string,uint64_t> refnametoid;
		std::vector<uint64_t> refoff;
		std::string refdata;
		loadfasta(reffa,refnametoid,refoff,refdata);

		std::map<std::string,uint64_t> readnametoid;
		std::vector<uint64_t> readoff;
		std::string readdata;
		loadfasta(rnfa,readnametoid,readoff,readdata);
		#endif

		#if 0
		{
			libmaus2::bambam::BamDecoder dec_a(fn_a);

			libmaus2::dazzler::align::AlignmentFileRegion::unique_ptr_type Palgn(libmaus2::dazzler::align::OverlapIndexer::openAlignmentFileWithoutIndex(ovl_a));

			libmaus2::dazzler::align::Overlap OVL;
			std::vector<libmaus2::dazzler::align::Overlap> VOVL;
			int64_t prevb = -1;
			int64_t rid = -1;
			int64_t rlen = 0;
			double gfrac = 0;
			double gcnt = 0;

			while ( Palgn->getNextOverlap(OVL) )
			{
				assert ( OVL.bread >= prevb );

				if ( OVL.bread != prevb )
				{
					if ( VOVL.size() )
					{
						gfrac += handleBlock(VOVL,dec_a,rid,rlen);
						gcnt += 1;
					}
					VOVL.resize(0);
				}

				VOVL.push_back(OVL);
				prevb = OVL.bread;
			}

			if ( VOVL.size() )
			{
				gfrac += handleBlock(VOVL,dec_a,rid,rlen);
				gcnt += 1;
			}
			VOVL.resize(0);

			std::cerr << "global " << gfrac / gcnt << std::endl;
		}
		#endif

		libmaus2::bambam::BamDecoder dec_a(fn_a);
		libmaus2::bambam::BamHeader const & header_a = dec_a.getHeader();
		libmaus2::bambam::BamDecoder dec_b(fn_b);
		// libmaus2::bambam::BamHeader const & header_b = dec_b.getHeader();

		struct ConditionalBamOutput
		{
			libmaus2::bambam::BamWriter::unique_ptr_type pout;

			ConditionalBamOutput(std::string const & fn, libmaus2::bambam::BamHeader const & rheader, bool const open)
			{
				if ( open )
				{
					libmaus2::bambam::BamWriter::unique_ptr_type tout(
						new libmaus2::bambam::BamWriter(fn,rheader)
					);
					pout = UNIQUE_PTR_MOVE(tout);
				}
			}

			void writeAlignment(libmaus2::bambam::BamAlignment const & algn)
			{
				if ( pout )
					pout->writeAlignment(algn);
			}
		};

		MissingOutput::unique_ptr_type nooverlapOut(new MissingOutput("no_overlap_out",header_a,writealgncat));
		ConditionalBamOutput primaryCrossedBW("primary_crossed.bam",header_a,writealgncat);
		// libmaus2::bambam::BamWriter::unique_ptr_type primaryCrossedBW(new libmaus2::bambam::BamWriter("primary_crossed.bam",header_a));
		ConditionalBamOutput anyCrossedBW("any_crossed.bam",header_a,writealgncat);
		//libmaus2::bambam::BamWriter::unique_ptr_type anyCrossedBW(new libmaus2::bambam::BamWriter("any_crossed.bam",header_a));
		ConditionalBamOutput primaryOverlapBW("primary_overlap.bam",header_a,writealgncat);
		//libmaus2::bambam::BamWriter::unique_ptr_type primaryOverlapBW(new libmaus2::bambam::BamWriter("primary_overlap.bam",header_a));
		ConditionalBamOutput anyOverlapBW("any_overlap.bam",header_a,writealgncat);
		// libmaus2::bambam::BamWriter::unique_ptr_type anyOverlapBW(new libmaus2::bambam::BamWriter("any_overlap.bam",header_a));
		ConditionalBamOutput primaryNoCrossBW("primary_no_cross.bam",header_a,writealgncat);
		// libmaus2::bambam::BamWriter::unique_ptr_type primaryNoCrossBW(new libmaus2::bambam::BamWriter("primary_no_cross.bam",header_a));
		ConditionalBamOutput primaryNoOverlapBW("primary_no_overlap.bam",header_a,writealgncat);
		// libmaus2::bambam::BamWriter::unique_ptr_type primaryNoOverlapBW(new libmaus2::bambam::BamWriter("primary_no_overlap.bam",header_a));

		libmaus2::lcs::NNPLocalAligner LA(6 /* bucket log */,14 /* k */,256*1024 /* max matches */,30 /* min band score */,50 /* min length */);

		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;

		BamPeeker apeeker(dec_a);
		BamPeeker bpeeker(dec_b);

		std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> Valgn_a;
		std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> Valgn_b;

		apeeker.getNextRange(Valgn_a);
		FAG.getNext();
		bpeeker.getNextRange(Valgn_b);

		uint64_t skipa = 0;
		uint64_t skipb = 0;

		uint64_t g_anycross = 0;
		uint64_t g_nocross = 0;
		uint64_t g_anyoverlap = 0;
		uint64_t g_nooverlap = 0;

		uint64_t g_primaryanycross = 0;
		uint64_t g_primarynocross = 0;
		uint64_t g_primaryanyoverlap = 0;
		uint64_t g_primarynooverlap = 0;

		uint64_t g_read_cross_bases = 0;
		uint64_t g_read_overlapbases = 0;
		uint64_t g_read_cross_bases_primary = 0;
		uint64_t g_read_overlapbases_primary = 0;
		uint64_t g_read_allbases = 0;

		uint64_t g_ref_cross_bases = 0;
		uint64_t g_ref_overlapbases = 0;
		uint64_t g_ref_cross_bases_primary = 0;
		uint64_t g_ref_overlapbases_primary = 0;
		uint64_t g_ref_allbases = 0;

		uint64_t unmapped = 0;

		libmaus2::bambam::StrCmpNum const strcmpnum;
		uint64_t c = 0;
		double mdiamacc = 0;
		double mdiamcnt = 0;

		double mdiamsum = 0;

		// uint64_t aid = 0;

		while ( Valgn_a.size() && Valgn_b.size() )
		{
			std::string const a_name = Valgn_a.front()->getName();
			std::string const b_name = Valgn_b.front()->getName();

			std::string shorter_name, longer_name;

			if ( a_name.size() <= b_name.size() )
			{
				shorter_name = a_name;
				longer_name = b_name;
			}
			else
			{
				shorter_name = b_name;
				longer_name = a_name;
			}

			std::string const longer_prefix = longer_name.substr(0,shorter_name.size());

			if ( longer_prefix == shorter_name )
			{
				assert ( Valgn_a.size() == 1 );
				libmaus2::bambam::BamAlignment const & a_algn = *(Valgn_a.front());

				libmaus2::geometry::RangeSet<RepeatLine> const * RS =
					(
						remapM.find(a_algn.getRefID()) != remapM.end()
					)
					? remapM.find(a_algn.getRefID())->second.get() : 0;

				if ( FAG.p )
				{
					assert ( FAG.pattern.getShortStringId() == a_algn.getName() );
					// std::cerr << a_algn.getName() << " " << FAG.pattern.sid << std::endl;
				}

				std::string readdata = a_algn.getRead();
				for ( uint64_t i = 0; i < readdata.size(); ++i )
					readdata[i] = libmaus2::fastx::mapChar(readdata[i]);

				uint64_t const numcig = libmaus2::bambam::BamAlignmentDecoderBase::getCigarOperations(Valgn_a[0]->D.begin(),cigop);

				libmaus2::lcs::AlignmentTraceContainer ATC;
				libmaus2::bambam::CigarStringParser::cigarToTrace(cigop.begin(),cigop.begin()+numcig,ATC);

				std::vector < std::pair<uint64_t,uint64_t> > kmatches = ATC.getKMatchOffsets(
					mapk,a_algn.getPos() - a_algn.getFrontDel(),a_algn.getFrontSoftClipping()/*offb*/
				);


				std::pair<uint64_t,uint64_t> P_PP = ATC.prefixPositive();
				std::pair<uint64_t,uint64_t> P_SP = ATC.suffixPositive();

				// error intervals in reference coordinates
				char const * c_eintv = a_algn.getAuxString("er");
				// error intervals in query coordinates
				//char const * c_rintv = a_algn.getAuxString("ee");

				std::vector<Intv> Vintv = c_eintv ? Intv::parse(std::string(c_eintv)) : std::vector<Intv>();

				#if 0
				// std::vector<Intv> Rintv = c_rintv ? Intv::parse(std::string(c_rintv)) : std::vector<Intv>();
				std::vector<Intv> Rintv = Intv::computeRIntv(Vintv, a_algn.D.begin());

				// shift error intervals in query coordinates
				for ( uint64_t i = 0; i < Rintv.size(); ++i )
				{
					Rintv[i].intv.from += a_algn.getFrontSoftClipping();
					Rintv[i].intv.to += a_algn.getFrontSoftClipping();
				}
				#endif

				uint64_t frontdel = 0;
				uint64_t zzz = 0;
				while (
					zzz < numcig
					&&
					cigop[zzz].first != libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH
					&&
					cigop[zzz].first != libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL
					&&
					cigop[zzz].first != libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF
				)
				{
					switch ( cigop[zzz].first )
					{
						case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
							frontdel += cigop[zzz].second;
							break;
					}

					zzz++;
				}

				// first reference base used in alignment
				int64_t const refstartpos = a_algn.getPos() - frontdel;

				// compute absolute intervals on reference
				for ( uint64_t i = 0; i < Vintv.size(); ++i )
				{
					Vintv[i].intv.from += refstartpos;
					Vintv[i].intv.to += refstartpos;
				}

				// perform self overlap analysis
				std::vector< libmaus2::math::IntegerInterval<int64_t> > TIV;
				if ( reference.size() )
				{
					assert ( a_algn.getRefID() >= 0 );
					assert ( Mref.find(a_algn.getRefID()) != Mref.end() );
					std::string const & ref = Mref.find(a_algn.getRefID())->second;
					libmaus2::math::IntegerInterval<int64_t> RI = a_algn.getReferenceInterval();
					RI.from -= frontdel;
					RI.to -= frontdel;
					assert ( RI.from >= 0 );
					assert ( RI.to <= static_cast<int64_t>(ref.size()) );

					std::string const sub = ref.substr(RI.from,RI.to+1-RI.from);

					std::vector< std::pair<libmaus2::lcs::NNPAlignResult,libmaus2::lcs::NNPTraceContainer::shared_ptr_type> > V = LA.align(
						sub.begin(),sub.end(),
						sub.begin(),sub.end()
					);

					// std::cerr << "V.size()=" << V.size() << std::endl;

					for ( uint64_t i = 0; i < V.size(); ++i )
						if ( V[i].first.bbpos < V[i].first.aepos )
						{
							TIV.push_back(libmaus2::math::IntegerInterval<int64_t>(RI.from + V[i].first.abpos,RI.from + V[i].first.bepos));
							// std::cerr << V[i].first << std::endl;
						}
					TIV = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(TIV);

					#if 0
					for ( uint64_t i = 0; i < TIV.size(); ++i )
						std::cerr << TIV[i] << std::endl;
					#endif

					LA.returnAlignments(V);
				}


				if ( FAG.p )
				{
					std::string & sid = FAG.pattern.sid;

					if ( FAG.pattern.sid.find("TAND") == std::string::npos )
					{
						std::ostringstream repaddstr;

						for ( uint64_t i = 0; i < TIV.size(); ++i )
							repaddstr << " TAND[" << TIV[i] << "]";

						sid += repaddstr.str();
					}
				}

				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_mapped;
				for ( uint64_t z = 0; z < Valgn_b.size(); ++z )
					if ( Valgn_b[z]->isMapped() )
						b_mapped.push_back(Valgn_b[z]);

				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_primary;
				for ( uint64_t z = 0; z < b_mapped.size(); ++z )
					if (
						!b_mapped[z]->isSecondary()
					)
						b_primary.push_back(b_mapped[z]);
				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_primary_correct_seq;
				for ( uint64_t z = 0; z < b_primary.size(); ++z )
					if (
						b_primary[z]->getRefID() == a_algn.getRefID()
					)
						b_primary_correct_seq.push_back(b_primary[z]);
				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_primary_correct_seq_correct_strand;
				for ( uint64_t z = 0; z < b_primary_correct_seq.size(); ++z )
					if (
						b_primary_correct_seq[z]->isReverse() == a_algn.isReverse()
					)
						b_primary_correct_seq_correct_strand.push_back(b_primary_correct_seq[z]);

				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vreadmapped;
				for ( uint64_t z = 0; z < b_mapped.size(); ++z )
					Vreadmapped.push_back(b_mapped[z]->getCoveredReadInterval());
				Vreadmapped = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vreadmapped);
				int64_t const slen = a_algn.getLseq() - ( a_algn.getFrontSoftClipping() + a_algn.getBackSoftClipping() );
				int64_t mdiam = 0;
				for ( uint64_t z = 0; z < Vreadmapped.size(); ++z )
					mdiam += Vreadmapped[z].diameter();
				double const mdiamfrac = static_cast<double>(mdiam) / slen;
				mdiamsum += mdiam;
				mdiamacc += mdiamfrac;
				mdiamcnt += 1;

				// std::cerr << "mdiamfrac=" << mdiamfrac << std::endl;

				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_correct_seq;
				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_incorrect_seq;
				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_correct_seq_correct_strand;
				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_correct_seq_incorrect_strand;

				if ( b_mapped.size() )
				{
					for ( uint64_t z = 0; z < b_mapped.size(); ++z )
						if ( b_mapped[z]->getRefID() == a_algn.getRefID() )
							b_correct_seq.push_back(b_mapped[z]);
						else
							b_incorrect_seq.push_back(b_mapped[z]);

					if ( b_correct_seq.size() )
					{
						for ( uint64_t z = 0; z < b_correct_seq.size(); ++z )
							if ( b_correct_seq[z]->isReverse() == a_algn.isReverse() )
								b_correct_seq_correct_strand.push_back(b_correct_seq[z]);
							else
								b_correct_seq_incorrect_strand.push_back(b_correct_seq[z]);
					}
				}
				else
				{
					unmapped += 1;
				}

				bool anycross = false;
				bool anyoverlap = false;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Voverlap_read;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vcross_read;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Voverlap_ref;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vcross_ref;

				for ( uint64_t z = 0; z < b_correct_seq_correct_strand.size(); ++z )
					if ( b_correct_seq_correct_strand[z]->getRead() != a_algn.getRead() )
					{
						std::string::size_type const p = a_algn.getRead().find(b_correct_seq_correct_strand[z]->getRead());

						if ( p != std::string::npos )
						{
							std::string cigstr = b_correct_seq_correct_strand[z]->getCigarString();
							std::ostringstream ostr;
							ostr << p << 'H';
							ostr << cigstr;
							ostr << a_algn.getRead().size() - b_correct_seq_correct_strand[z]->getRead().size() - p << 'H';
							b_correct_seq_correct_strand[z]->replaceCigarString(ostr.str());
						}
					}

				for ( uint64_t z = 0; z < b_correct_seq_correct_strand.size(); ++z )
				{
					#if 0
					std::cerr
						<< a_algn.getName() << " "
						<< b_correct_seq_correct_strand[z]->getName() << " "
						<< a_algn.getReferenceInterval() << " " << b_correct_seq_correct_strand[z]->getReferenceInterval() << std::endl;
					#endif

					// if we have any overlap between the correct reference region and the mapped one
					if ( !
						a_algn.getReferenceInterval().intersection(
							b_correct_seq_correct_strand[z]->getReferenceInterval()
						).isEmpty()
					)
					{
						anyoverlap = true;

						// crossing between real alignment and computed one?
						bool const cross = libmaus2::bambam::BamAlignment::cross(a_algn,*(b_correct_seq_correct_strand[z]));
						anycross = anycross || cross;

						anyOverlapBW.writeAlignment(*b_correct_seq_correct_strand[z]);
						if ( cross )
							anyCrossedBW.writeAlignment(*b_correct_seq_correct_strand[z]);

						if ( cross )
							Vcross_ref.push_back(b_correct_seq_correct_strand[z]->getReferenceInterval());
						Voverlap_ref.push_back(b_correct_seq_correct_strand[z]->getReferenceInterval());

						if ( cross )
							Vcross_read.push_back(b_correct_seq_correct_strand[z]->getCoveredReadInterval());
						Voverlap_read.push_back(b_correct_seq_correct_strand[z]->getCoveredReadInterval());
					}
				}

				bool primaryanyoverlap = false;
				bool primaryanycross = false;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vprimaryoverlap_read;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vprimarycross_read;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vprimaryoverlap_ref;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vprimarycross_ref;

				for ( uint64_t z = 0; z < b_primary_correct_seq_correct_strand.size(); ++z )
				{
					if ( !
						a_algn.getReferenceInterval().intersection(
							b_primary_correct_seq_correct_strand[z]->getReferenceInterval()
						).isEmpty()
					)
					{
						primaryanyoverlap = true;

						bool const cross = libmaus2::bambam::BamAlignment::cross(a_algn,*(b_primary_correct_seq_correct_strand[z]));
						primaryanycross = primaryanycross || cross;

						primaryOverlapBW.writeAlignment(*b_primary_correct_seq_correct_strand[z]);
						if ( cross )
							primaryCrossedBW.writeAlignment(*b_primary_correct_seq_correct_strand[z]);

						if ( cross )
							Vprimarycross_ref.push_back(b_primary_correct_seq_correct_strand[z]->getReferenceInterval());
						Vprimaryoverlap_ref.push_back(b_primary_correct_seq_correct_strand[z]->getReferenceInterval());

						if ( cross )
							Vprimarycross_read.push_back(b_primary_correct_seq_correct_strand[z]->getCoveredReadInterval());
						Vprimaryoverlap_read.push_back(b_primary_correct_seq_correct_strand[z]->getCoveredReadInterval());

					}
				}

				if ( ! primaryanyoverlap )
					primaryNoOverlapBW.writeAlignment(a_algn);
				if ( ! primaryanycross )
				{
					primaryNoCrossBW.writeAlignment(a_algn);

					if ( RS )
					{
						uint64_t const left = a_algn.getPos();
						uint64_t const right = left + a_algn.getReferenceLength();
						std::cerr << "MISSED PRIMARY [" << left << "," << right << ") ";

						double repcnt = 0;
						double repavgcnt = 0;
						uint64_t repavgden = 0;


						RepeatLine RL(a_algn.getRefID(),left,right-left,0.0/*e */);

						std::vector<RepeatLine const *> const VV = RS->search(RL);
						// libmaus2::geometry::RangeSet<RepeatLine> const * RS

						std::cerr << VV.size() << " ";

						for ( uint64_t i = 0; i < VV.size(); ++i )
						{
							libmaus2::math::IntegerInterval<int64_t> O = VV[i]->getOverlap(a_algn);
							int64_t const diam = O.diameter();
							repcnt += diam;
							repavgcnt += VV[i]->e * diam;
							repavgden += diam;

							std::cerr << (*(VV[i]));
						}

						std::cerr << " rep=" << repcnt / (right-left) << " e=" << repavgcnt/repavgden << std::endl;
					}

				}

				if ( primaryanyoverlap )
					g_primaryanyoverlap += 1;
				else
					g_primarynooverlap += 1;
				if ( primaryanycross )
					g_primaryanycross += 1;
				else
					g_primarynocross += 1;

				// merge intervals on read (overlap + cross)
				Voverlap_read = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Voverlap_read);
				Vcross_read = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vcross_read);

				// merge intervals on reference (overlap + cross)
				Voverlap_ref = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Voverlap_ref);
				Vcross_ref = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vcross_ref);

				// merge intervals on read (overlap + cross)
				Vprimaryoverlap_read = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vprimaryoverlap_read);
				Vprimarycross_read = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vprimarycross_read);

				// merge intervals on reference (overlap + cross)
				Vprimaryoverlap_ref = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vprimaryoverlap_ref);
				Vprimarycross_ref = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vprimarycross_ref);

				std::vector< libmaus2::math::IntegerInterval<int64_t> >
					notCrossed_read = libmaus2::math::IntegerInterval<int64_t>::difference(a_algn.getCoveredReadInterval(),Vcross_read);
				std::vector< libmaus2::math::IntegerInterval<int64_t> >
					notCrossed_ref = libmaus2::math::IntegerInterval<int64_t>::difference(a_algn.getReferenceInterval(),Vcross_ref);

				std::vector< libmaus2::math::IntegerInterval<int64_t> >
					primarynotCrossed_read = libmaus2::math::IntegerInterval<int64_t>::difference(a_algn.getCoveredReadInterval(),Vprimarycross_read);
				std::vector< libmaus2::math::IntegerInterval<int64_t> >
					primarynotCrossed_ref = libmaus2::math::IntegerInterval<int64_t>::difference(a_algn.getReferenceInterval(),Vprimarycross_ref);

				if ( notCrossed_ref.size() )
				{
					#if 0
					for ( uint64_t i = 0; i < Vintv.size(); ++i )
						std::cerr << Vintv[i] << std::endl;
					#endif
					libmaus2::math::IntegerInterval<int64_t> const RI = a_algn.getReferenceInterval();

					for ( uint64_t j = 0; j < notCrossed_ref.size(); ++j )
					{
						bool const front = (RI.from == notCrossed_ref[j].from);
						bool const back = (RI.to == notCrossed_ref[j].to);

						// allow front/back clipping of this number of bases without reporting them as missing
						uint64_t const frontbackclipallow = 20;
						// allow this many low error bases to be missed without reporting
						//uint64_t const lowerrallow = 150;

						if ( (front || back) && (notCrossed_ref[j].diameter() <= static_cast<int64_t>(frontbackclipallow)) )
							continue;

						// get minimum error rate over all relevant intervals
						double minerr = std::numeric_limits<double>::max();
						for ( uint64_t i = 0; i < Vintv.size(); ++i )
							if ( ! notCrossed_ref[j].intersection(Vintv[i].intv).isEmpty() )
								if ( Vintv[i].erate < minerr )
									minerr = Vintv[i].erate;

						// if all error rates are too high
						if ( minerr >= errthres )
							continue;

						// count number of low error rate bases
						uint64_t lowerrbases = 0;
						for ( uint64_t i = 0; i < Vintv.size(); ++i )
							if ( ! notCrossed_ref[j].intersection(Vintv[i].intv).isEmpty() )
								if ( Vintv[i].erate <= errthres )
									lowerrbases += notCrossed_ref[j].intersection(Vintv[i].intv).diameter();

						// ignore if number is sufficiently low
						if ( lowerrbases <= lowerrallow )
							continue;

						std::cerr << "not crossed on reference " << notCrossed_ref[j] << " " << RI << " " << minerr << " " << P_PP.first << "," << P_SP.first << " lowerr=" << lowerrbases;

						if ( front )
							std::cerr << " front";
						if ( back )
							std::cerr << " back";

						for ( uint64_t i = 0; i < Vintv.size(); ++i )
							if ( ! notCrossed_ref[j].intersection(Vintv[i].intv).isEmpty() )
								std::cerr << " " << Vintv[i];

						std::cerr << std::endl;

						for ( uint64_t i = 0; i < TIV.size(); ++i )
							if ( ! TIV[i].intersection(notCrossed_ref[j]).isEmpty() )
								std::cerr << TIV[i] << std::endl;

						#if 0
						for ( uint64_t z = 0; z < b_mapped.size(); ++z )
						{
							::libmaus2::bambam::BamFormatAuxiliary aux;
							b_mapped[z]->formatAlignment(std::cerr,header_b,aux);
							std::cerr << std::endl;
						}
						#endif

						#if 0
						std::cerr << "crossed: ";
						for ( uint64_t i = 0; i < Vcross_ref.size(); ++i )
							std::cerr << Vcross_ref[i];
						std::cerr << std::endl;
						std::cerr << "overlap: ";
						for ( uint64_t i = 0; i < Vcross_ref.size(); ++i )
							std::cerr << Voverlap_ref[i];
						std::cerr << std::endl;
						#endif
					}
				}

				#if 0
				for ( uint64_t i = 0; i < notCrossed_read.size(); ++i )
				{
					std::cerr << "not crossed " << notCrossed_read[i] << std::endl;

					for ( uint64_t j = 0; j < Rintv.size(); ++j )
						if ( !notCrossed_read[i].intersection(Rintv[j].intv).isEmpty() )
							std::cerr << Rintv[j] << std::endl;
				}
				#endif

				for ( uint64_t z = 0; z < Voverlap_read.size(); ++z )
					Voverlap_read[z] = a_algn.getCoveredReadInterval().intersection(Voverlap_read[z]);
				for ( uint64_t z = 0; z < Vcross_read.size(); ++z )
					Vcross_read[z] = a_algn.getCoveredReadInterval().intersection(Vcross_read[z]);

				for ( uint64_t z = 0; z < Voverlap_ref.size(); ++z )
					Voverlap_ref[z] = a_algn.getReferenceInterval().intersection(Voverlap_ref[z]);
				for ( uint64_t z = 0; z < Vcross_ref.size(); ++z )
					Vcross_ref[z] = a_algn.getReferenceInterval().intersection(Vcross_ref[z]);

				for ( uint64_t z = 0; z < Vprimaryoverlap_read.size(); ++z )
					Vprimaryoverlap_read[z] = a_algn.getCoveredReadInterval().intersection(Vprimaryoverlap_read[z]);
				for ( uint64_t z = 0; z < Vprimarycross_read.size(); ++z )
					Vprimarycross_read[z] = a_algn.getCoveredReadInterval().intersection(Vprimarycross_read[z]);

				for ( uint64_t z = 0; z < Vprimaryoverlap_ref.size(); ++z )
					Vprimaryoverlap_ref[z] = a_algn.getReferenceInterval().intersection(Vprimaryoverlap_ref[z]);
				for ( uint64_t z = 0; z < Vprimarycross_ref.size(); ++z )
					Vprimarycross_ref[z] = a_algn.getReferenceInterval().intersection(Vprimarycross_ref[z]);

				uint64_t ovl_read_bases = 0;
				for ( uint64_t z = 0; z < Voverlap_read.size(); ++z )
					ovl_read_bases += Voverlap_read[z].diameter();
				uint64_t read_cross_bases = 0;
				for ( uint64_t z = 0; z < Vcross_read.size(); ++z )
					read_cross_bases += Vcross_read[z].diameter();

				uint64_t ovl_ref_bases = 0;
				for ( uint64_t z = 0; z < Voverlap_ref.size(); ++z )
					ovl_ref_bases += Voverlap_ref[z].diameter();
				uint64_t ref_cross_bases = 0;
				for ( uint64_t z = 0; z < Vcross_ref.size(); ++z )
					ref_cross_bases += Vcross_ref[z].diameter();

				uint64_t ovl_read_bases_primary = 0;
				for ( uint64_t z = 0; z < Vprimaryoverlap_read.size(); ++z )
					ovl_read_bases_primary += Vprimaryoverlap_read[z].diameter();
				uint64_t read_cross_bases_primary = 0;
				for ( uint64_t z = 0; z < Vprimarycross_read.size(); ++z )
					read_cross_bases_primary += Vprimarycross_read[z].diameter();

				uint64_t ovl_ref_bases_primary = 0;
				for ( uint64_t z = 0; z < Vprimaryoverlap_ref.size(); ++z )
					ovl_ref_bases_primary += Vprimaryoverlap_ref[z].diameter();
				uint64_t ref_cross_bases_primary = 0;
				for ( uint64_t z = 0; z < Vprimarycross_ref.size(); ++z )
					ref_cross_bases_primary += Vprimarycross_ref[z].diameter();

				if ( anycross )
					g_anycross += 1;
				else
				{
					g_nocross += 1;

					if ( RS )
					{
						uint64_t const left = a_algn.getPos();
						uint64_t const right = left + a_algn.getReferenceLength();
						std::cerr << "MISSED ANY [" << left << "," << right << ") ";

						double repcnt = 0;
						double repavgcnt = 0;
						uint64_t repavgden = 0;

						RepeatLine RL(a_algn.getRefID(),left,right-left,0.0/*e */);

						std::vector<RepeatLine const *> const VV = RS->search(RL);
						// libmaus2::geometry::RangeSet<RepeatLine> const * RS

						std::cerr << VV.size() << " ";

						for ( uint64_t i = 0; i < VV.size(); ++i )
						{
							libmaus2::math::IntegerInterval<int64_t> O = VV[i]->getOverlap(a_algn);
							int64_t const diam = O.diameter();
							repcnt += diam;
							repavgcnt += VV[i]->e * diam;
							repavgden += diam;

							std::cerr << (*(VV[i]));
						}

						std::cerr << " rep=" << repcnt / (right-left) << " e=" << repavgcnt/repavgden << std::endl;

						if ( Prank )
						{
							std::cerr << "match ";
							for ( uint64_t i = 0; i < kmatches.size(); ++i )
							{
								uint64_t const rpos = kmatches[i].second;
								std::pair<uint64_t,uint64_t> const P = Prank->backwardSearch(readdata.begin()+rpos,mapk);
								std::cerr << P.second-P.first << ((i+1 < kmatches.size())?",":"");
							}
							std::cerr << std::endl;
						}
					}
				}

				if ( anyoverlap )
					g_anyoverlap += 1;
				else
					g_nooverlap += 1;

				if ( ! anyoverlap )
				{
					if ( FAG.p )
						nooverlapOut->put(a_algn,FAG.pattern);
					else
						nooverlapOut->put(a_algn);
				}

				g_read_cross_bases += read_cross_bases;
				g_read_overlapbases += ovl_read_bases;

				g_ref_cross_bases += ref_cross_bases;
				g_ref_overlapbases += ovl_ref_bases;

				g_read_cross_bases_primary += read_cross_bases_primary;
				g_read_overlapbases_primary += ovl_read_bases_primary;

				g_ref_cross_bases_primary += ref_cross_bases_primary;
				g_ref_overlapbases_primary += ovl_ref_bases_primary;

				g_ref_allbases += a_algn.getReferenceLength();
				g_read_allbases += a_algn.getCoveredReadInterval().diameter();

				if ( ++c % verbmod == 0 )
				{
					std::cerr << "[D] "
						<< c << " "
						<< a_algn.getName()
						<< ","
						<< a_algn.isReverse()
						<< " " << a_algn.getReferenceLength() << " " << a_algn.getErrorRate()
						<< " correct seq/strand " << b_correct_seq_correct_strand.size()
						<< " wrong strand " << b_correct_seq_incorrect_strand.size()
						<< " wrong seq " << b_incorrect_seq.size()
						<< " cross " << anycross
						<< " overlap " << anyoverlap
						<< " pcross " << primaryanycross
						<< " poverlap " << primaryanyoverlap
						<< " read_cross_bases " << read_cross_bases << " " << (static_cast<double>(read_cross_bases) / a_algn.getCoveredReadInterval().diameter())
						<< " ovl_read_bases " << ovl_read_bases << " " << (static_cast<double>(ovl_read_bases) / a_algn.getCoveredReadInterval().diameter())
						<< " pcscs " << b_primary_correct_seq_correct_strand.size()
						<< " " << static_cast<double>(g_read_overlapbases) / g_read_allbases
						<< " " << static_cast<double>(g_read_cross_bases) / g_read_allbases
						<< " " << static_cast<double>(g_ref_overlapbases) / g_ref_allbases
						<< " " << static_cast<double>(g_ref_cross_bases) / g_ref_allbases
						<< " " << static_cast<double>(g_read_overlapbases_primary) / g_read_allbases
						<< " " << static_cast<double>(g_read_cross_bases_primary) / g_read_allbases
						<< " " << static_cast<double>(g_ref_overlapbases_primary) / g_ref_allbases
						<< " " << static_cast<double>(g_ref_cross_bases_primary) / g_ref_allbases
						<< " " << skipa
						<< " " << skipb
						<< " " << unmapped
						<< " " << g_anycross << " (" << static_cast<double>(g_anycross)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< " " << g_nocross << " (" << static_cast<double>(g_nocross)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< " " << g_anyoverlap << " (" << static_cast<double>(g_anyoverlap)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< " " << g_nooverlap << " (" << static_cast<double>(g_nooverlap)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< " " << g_primaryanycross << " (" << static_cast<double>(g_primaryanycross)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< " " << g_primarynocross << " (" << static_cast<double>(g_primarynocross)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< " " << g_primaryanyoverlap << " (" << static_cast<double>(g_primaryanyoverlap)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< " " << g_primarynooverlap << " (" << static_cast<double>(g_primarynooverlap)/static_cast<double>(g_anycross+g_nocross+skipa+skipb+unmapped) << ")"
						<< std::endl;
				}

				apeeker.getNextRange(Valgn_a);
				FAG.getNext();
				bpeeker.getNextRange(Valgn_b);
			}
			else
			{
				int const cmp = strcmpnum.strcmpnum(Valgn_a.front()->getName(),Valgn_b.front()->getName());

				assert ( cmp != 0 );

				if ( cmp < 0 )
				{
					//std::cerr << "skip a " << Valgn_a[0]->getErrorRate() << " " << Valgn_a[0]->getReferenceLength() << std::endl;
					g_read_allbases += Valgn_a[0]->getCoveredReadInterval().diameter();
					g_ref_allbases += Valgn_a[0]->getReferenceLength();
					mdiamcnt += 1;
					apeeker.getNextRange(Valgn_a);
					FAG.getNext();
					++skipa;
				}
				else if ( cmp > 0 )
				{
					//std::cerr << "skip b" << std::endl;
					bpeeker.getNextRange(Valgn_b);
					++skipb;
				}
			}
		}

		std::cerr << "[S] skip_a=" << skipa << std::endl;
		std::cerr << "[S] skip_b=" << skipb << std::endl;
		std::cerr << "[S] unmapped=" << unmapped << std::endl;
		std::cerr << "[S] g_anycross=" << g_anycross << std::endl;
		std::cerr << "[S] g_nocross=" << g_nocross << std::endl;
		std::cerr << "[S] g_anyoverlap=" << g_anyoverlap << std::endl;
		std::cerr << "[S] g_nooverlap=" << g_nooverlap << std::endl;
		std::cerr << "[S] g_primaryanycross=" << g_primaryanycross << std::endl;
		std::cerr << "[S] g_primarynocross=" << g_primarynocross << std::endl;
		std::cerr << "[S] g_primaryanyoverlap=" << g_primaryanyoverlap << std::endl;
		std::cerr << "[S] g_primarynooverlap=" << g_primarynooverlap << std::endl;
		if ( g_read_allbases )
		{
			std::cerr << "[S] read overlap fraction=" << static_cast<double>(g_read_overlapbases) / g_read_allbases << std::endl;
			std::cerr << "[S] read cross fraction=" << static_cast<double>(g_read_cross_bases) / g_read_allbases << std::endl;
		}
		if ( g_ref_allbases )
		{
			std::cerr << "[S] ref overlap fraction=" << static_cast<double>(g_ref_overlapbases) / g_ref_allbases << std::endl;
			std::cerr << "[S] ref cross fraction=" << static_cast<double>(g_ref_cross_bases) / g_ref_allbases << std::endl;
		}
		if ( g_read_allbases )
		{
			std::cerr << "[S] primary read overlap fraction=" << static_cast<double>(g_read_overlapbases_primary) / g_read_allbases << std::endl;
			std::cerr << "[S] primary read cross fraction=" << static_cast<double>(g_read_cross_bases_primary) / g_read_allbases << std::endl;
		}
		if ( g_ref_allbases )
		{
			std::cerr << "[S] primary ref overlap fraction=" << static_cast<double>(g_ref_overlapbases_primary) / g_ref_allbases << std::endl;
			std::cerr << "[S] primary ref cross fraction=" << static_cast<double>(g_ref_cross_bases_primary) / g_ref_allbases << std::endl;
		}

		std::cerr << "[S] mdiam avg " << mdiamacc/mdiamcnt << std::endl;
		std::cerr << "[S] mdiam sum avg " << mdiamsum / g_read_allbases << std::endl;

		std::string nooverlapfastafn = nooverlapOut->fastafilename;
		nooverlapOut.reset();
		if ( ! (FAG.p) )
		{
			libmaus2::aio::FileRemoval::removeFile(nooverlapfastafn);
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
