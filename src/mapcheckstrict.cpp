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

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgParser const arg(argc,argv);
		std::string const fn_a = arg[0];
		std::string const fn_b = arg[1];
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
		libmaus2::bambam::BamHeader const & header_b = dec_b.getHeader();

		BamPeeker apeeker(dec_a);
		BamPeeker bpeeker(dec_b);

		std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> Valgn_a;
		std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> Valgn_b;

		apeeker.getNextRange(Valgn_a);
		bpeeker.getNextRange(Valgn_b);

		uint64_t skipa = 0;
		uint64_t skipb = 0;
		uint64_t g_anycross = 0;
		uint64_t g_nocross = 0;
		uint64_t g_anyoverlap = 0;
		uint64_t g_nooverlap = 0;

		uint64_t g_read_cross_bases = 0;
		uint64_t g_read_overlapbases = 0;
		uint64_t g_read_allbases = 0;

		uint64_t g_ref_cross_bases = 0;
		uint64_t g_ref_overlapbases = 0;
		uint64_t g_ref_allbases = 0;

		libmaus2::bambam::StrCmpNum const strcmpnum;
		uint64_t c = 0;
		uint64_t verbmod = 1024;
		double mdiamacc = 0;
		double mdiamcnt = 0;

		double mdiamsum = 0;

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

				char const * c_eintv = a_algn.getAuxString("er");
				char const * c_rintv = a_algn.getAuxString("ee");

				std::vector<Intv> Vintv = c_eintv ? Intv::parse(std::string(c_eintv)) : std::vector<Intv>();
				// std::vector<Intv> Rintv = c_rintv ? Intv::parse(std::string(c_rintv)) : std::vector<Intv>();
				std::vector<Intv> Rintv = Intv::computeRIntv(Vintv, a_algn.D.begin());

				for ( uint64_t i = 0; i < Rintv.size(); ++i )
				{
					Rintv[i].intv.from += a_algn.getFrontSoftClipping();
					Rintv[i].intv.to += a_algn.getFrontSoftClipping();
				}

				#if 0
				uint64_t frontdel = 0;
				while ( zz < numcig && cigop[zz].first == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL )
					frontdel += cigop[zz++].second;

				int64_t const intvshift = a_algn.getPos() - frontdel;

				for ( uint64_t i = 0; i < Vintv.size(); ++i )
				{
					Vintv[i].intv.from += intvshift;
					Vintv[i].intv.to += intvshift;
				}
				#endif

				std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type> b_mapped;
				for ( uint64_t z = 0; z < Valgn_b.size(); ++z )
					if ( Valgn_b[z]->isMapped() )
						b_mapped.push_back(Valgn_b[z]);

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

				bool anycross = false;
				bool anyoverlap = false;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Voverlap_read;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vcross_read;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Voverlap_ref;
				std::vector< libmaus2::math::IntegerInterval<int64_t> > Vcross_ref;

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

						if ( cross )
							Vcross_ref.push_back(b_correct_seq_correct_strand[z]->getReferenceInterval());
						Voverlap_ref.push_back(b_correct_seq_correct_strand[z]->getReferenceInterval());

						if ( cross )
							Vcross_read.push_back(b_correct_seq_correct_strand[z]->getCoveredReadInterval());
						Voverlap_read.push_back(b_correct_seq_correct_strand[z]->getCoveredReadInterval());
					}
				}

				Voverlap_read = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Voverlap_read);
				Vcross_read = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vcross_read);

				Voverlap_ref = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Voverlap_ref);
				Vcross_ref = libmaus2::math::IntegerInterval<int64_t>::mergeTouchingOrOverlapping(Vcross_ref);

				std::vector< libmaus2::math::IntegerInterval<int64_t> >
					notCrossed_read = libmaus2::math::IntegerInterval<int64_t>::difference(
						a_algn.getCoveredReadInterval(),
						Vcross_read
					);

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

				if ( anycross )
					g_anycross += 1;
				else
					g_nocross += 1;

				if ( anyoverlap )
					g_anyoverlap += 1;
				else
					g_nooverlap += 1;

				g_read_cross_bases += read_cross_bases;
				g_read_overlapbases += ovl_read_bases;

				g_ref_cross_bases += ref_cross_bases;
				g_ref_overlapbases += ovl_ref_bases;

				g_ref_allbases += a_algn.getReferenceLength();
				g_read_allbases += a_algn.getCoveredReadInterval().diameter();

				if ( ++c % verbmod == 0 )
				{
					std::cerr << "[D] "
						<< c << " "
						<< a_algn.getName() << " " << a_algn.getReferenceLength() << " " << a_algn.getErrorRate()
						<< " correct seq/strand " << b_correct_seq_correct_strand.size()
						<< " wrong strand " << b_correct_seq_incorrect_strand.size()
						<< " wrong seq " << b_incorrect_seq.size()
						<< " cross " << anycross
						<< " overlap " << anyoverlap
						<< " read_cross_bases " << read_cross_bases << " " << (static_cast<double>(read_cross_bases) / a_algn.getCoveredReadInterval().diameter())
						<< " ovl_read_bases " << ovl_read_bases << " " << (static_cast<double>(ovl_read_bases) / a_algn.getCoveredReadInterval().diameter())
						<< " " << static_cast<double>(g_read_overlapbases) / g_read_allbases
						<< " " << static_cast<double>(g_read_cross_bases) / g_read_allbases
						<< " " << static_cast<double>(g_ref_overlapbases) / g_ref_allbases
						<< " " << static_cast<double>(g_ref_cross_bases) / g_ref_allbases
						<< std::endl;
				}

				apeeker.getNextRange(Valgn_a);
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
		std::cerr << "[S] g_anycross=" << g_anycross << std::endl;
		std::cerr << "[S] g_nocross=" << g_nocross << std::endl;
		std::cerr << "[S] g_anyoverlap=" << g_anyoverlap << std::endl;
		std::cerr << "[S] g_nooverlap=" << g_nooverlap << std::endl;
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

		std::cerr << "[S] mdiam avg " << mdiamacc/mdiamcnt << std::endl;
		std::cerr << "[S] mdiam sum avg " << mdiamsum / g_read_allbases << std::endl;

		#if 0
		std::cerr << "[V] Reading alignments...";
		while ( dec_a.readAlignment() )
			Valgn_a.push_back(dec_a.getAlignment().sclone());
		while ( dec_b.readAlignment() )
			Valgn_b.push_back(dec_b.getAlignment().sclone());
		std::cerr << "done." << std::endl;

		libmaus2::bambam::StrCmpNum const strcmpnum;
		::libmaus2::bambam::BamFormatAuxiliary aux;
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Acigop;
		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Bcigop;

		libmaus2::lcs::AlignmentTraceContainer ATC_a;
		libmaus2::lcs::AlignmentTraceContainer ATC_b;

		typedef std::vector<libmaus2::bambam::BamAlignment::shared_ptr_type>::const_iterator iterator;
		iterator a_itc = Valgn_a.begin();
		iterator b_itc = Valgn_b.begin();

		uint64_t numcross = 0;
		uint64_t numwrongstrand = 0;
		uint64_t numwrongrefid = 0;
		uint64_t numnocross = 0;
		uint64_t numnoprima = 0;
		uint64_t numnoprimb = 0;
		uint64_t skipa = 0;
		uint64_t skipb = 0;
		double fracsum = 0;

		while ( a_itc != Valgn_a.end() && b_itc != Valgn_b.end() )
		{
			std::string const a_name = (*a_itc)->getName();
			std::string const b_name = (*b_itc)->getName();
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

			// same name
			if ( longer_prefix == shorter_name )
			{
				std::cerr << std::string(80,'-') << std::endl;

				assert ( readnametoid.find(b_name) != readnametoid.end() );
				OVL.bread = readnametoid.find(b_name)->second;
				typedef std::vector<libmaus2::dazzler::align::Overlap>::const_iterator o_it;
				std::pair<o_it,o_it> const OP = std::equal_range(VO.begin(),VO.end(),OVL,OverlapBComparator());

				uint64_t const b_padreadlen = (readoff.at(OVL.bread+1)-readoff.at(OVL.bread))/2;
				uint64_t const b_readlen = b_padreadlen - 2;

				#if 0
				uint64_t cc = 0;
				for ( uint64_t i = 0; i < VO.size(); ++i )
					if ( VO[i].bread == OVL.bread )
						++cc;
				#endif

				iterator a_top = a_itc;
				while ( a_top != Valgn_a.end() && strcmpnum.strcmpnum((*a_itc)->getName(),(*a_top)->getName()) == 0 )
					++a_top;
				iterator b_top = b_itc;
				while ( b_top != Valgn_b.end() && strcmpnum.strcmpnum((*b_itc)->getName(),(*b_top)->getName()) == 0 )
					++b_top;

				std::cerr << (*a_itc)->getName() << " " << (a_top-a_itc) << " " << (b_top-b_itc) << " " << OP.second-OP.first << " " << OVL.bread << std::endl;

				//assert ( OP.second-OP.first == b_top-b_itc );

				for ( o_it it = OP.first; it != OP.second; ++it )
				{
					libmaus2::dazzler::align::Overlap const & OVL = *it;
					char const * adata = refdata.c_str() + refoff.at(OVL.aread) + 1;

					bool const isinv = OVL.isInverse();
					uint64_t const padreadlen = (readoff.at(OVL.bread+1)-readoff.at(OVL.bread))/2;
					uint64_t const readlen = padreadlen - 2;
					char const * bdata = readdata.c_str() + readoff.at(OVL.bread) + (isinv ? padreadlen+1 : 1);

					libmaus2::lcs::AlignmentTraceContainer ATC;
					libmaus2::lcs::NP aligner;

					libmaus2::dazzler::align::Overlap::computeTracePreMapped(
						OVL.path.path.begin(),OVL.path.path.size(),OVL.path.abpos,OVL.path.aepos,OVL.path.bbpos,OVL.path.bepos,
						reinterpret_cast<uint8_t const *>(adata),reinterpret_cast<uint8_t const *>(bdata),tspace,ATC,aligner);

					libmaus2::autoarray::AutoArray< std::pair<libmaus2::lcs::AlignmentTraceContainer::step_type,uint64_t> > Aopblocks;
					libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> Aop;

					uint64_t const numcig = libmaus2::bambam::CigarStringParser::traceToCigar(ATC,Aopblocks,Aop,0 /* hard clip left */,
						OVL.path.bbpos,readlen - OVL.path.bepos,0 /* hard clip right */
					);
					std::string const cigstr = libmaus2::bambam::CigarStringParser::cigarBlocksToString(Aop.begin(),numcig);
					// std::cerr << "OVL " << isinv << " "  << cigstr << std::endl;
				}

				#if 0
				for ( iterator it = b_itc; it != b_top; ++it )
				{
					std::cerr << (*it)->getCigarString() << std::endl;
				}
				#endif

				if ( a_top - a_itc != 1 )
				{
					std::cerr << "[E] ambiguous for " << a_name << " " << b_name << std::endl;
					assert ( a_top-a_itc == 1 );
				}

				std::vector<std::pair<int64_t,int64_t> > Amappos;
				(*a_itc)->getMappingPositionPairs(Amappos);
				std::sort(Amappos.begin(),Amappos.end());

				std::vector<std::pair<int64_t,int64_t> > ABmappos;

				bool anycross = false;
				uint64_t numcorrectseq = 0;
				int64_t mindist = std::numeric_limits<int64_t>::max();
				iterator best;
				int64_t bestdist = std::numeric_limits<int64_t>::max();

				for ( iterator it = b_itc; it != b_top; ++it )
				{
					libmaus2::bambam::BamAlignment::shared_ptr_type const & prim_a = *a_itc;
					libmaus2::bambam::BamAlignment::shared_ptr_type const & prim_b = *it;

					std::vector<std::pair<int64_t,int64_t> > Bmappos;

					bool const ccross = ::libmaus2::bambam::BamAlignment::cross(*prim_a,*prim_b);

					anycross = anycross || ccross;

					if ( ccross )
					{
						best = it;
						bestdist = -1;
						mindist = 0;

						prim_b->getMappingPositionPairs(Bmappos);
						std::sort(Bmappos.begin(),Bmappos.end());
						std::vector<std::pair<int64_t,int64_t> > inters;
						std::set_intersection(
							Amappos.begin(),Amappos.end(),
							Bmappos.begin(),Bmappos.end(),
							std::back_insert_iterator< std::vector<std::pair<int64_t,int64_t> > >(inters)
						);
						std::copy(inters.begin(),inters.end(),std::back_insert_iterator< std::vector<std::pair<int64_t,int64_t> > >(ABmappos));
						ABmappos.resize( std::unique(ABmappos.begin(),ABmappos.end()) - ABmappos.begin() );
					}
					else
					{
						if ( prim_a->getRefID() != prim_b->getRefID() )
						{

						}
						else if ( prim_a->isReverse() != prim_b->isReverse() )
						{

						}
						else
						{
							numcorrectseq += 1;

							libmaus2::math::IntegerInterval<int64_t> const aintv = prim_a->getReferenceInterval();
							libmaus2::math::IntegerInterval<int64_t> const bintv = prim_b->getReferenceInterval();

							if ( ! aintv.intersection(bintv).isEmpty() )
							{
								if ( bestdist > 0 )
								{
									mindist = 0;
									bestdist = 0;
									best = it;
								}
							}
							else
							{
								int64_t avdist;

								if ( aintv.from > bintv.to )
									avdist = aintv.from- bintv.to;
								else
									avdist = bintv.from - aintv.to;

								if ( avdist < mindist )
								{
									mindist = avdist;
									bestdist = avdist;
									best = it;
								}
							}

						}
					}

					#if 0
					std::cerr << "ccross=" << ccross << " " << prim_b->getCigarString() << std::endl;

					if ( prim_a->getRefID() != prim_b->getRefID() )
					{

					}
					else if ( prim_a->isReverse() != prim_b->isReverse() )
					{

					}
					else
					{
						numcorrectseq += 1;

						uint64_t const start_a_r = prim_a->getPos();
						uint64_t const start_a_q = keep_soft_clipping ? prim_a->getFrontSoftClipping() : 0;
						uint64_t const start_b_r = prim_b->getPos();
						uint64_t const start_b_q = keep_soft_clipping ? prim_b->getFrontSoftClipping() : 0;

						uint32_t const num_cig_a = prim_a->getCigarOperations(Acigop);
						uint32_t const num_cig_b = prim_b->getCigarOperations(Bcigop);

						libmaus2::bambam::CigarStringParser::cigarToTrace(Acigop.begin(),Acigop.begin()+num_cig_a,ATC_a,true /* ignore unknown */);
						libmaus2::bambam::CigarStringParser::cigarToTrace(Bcigop.begin(),Bcigop.begin()+num_cig_b,ATC_b,true /* ignore unknown */);

						uint64_t offset_a = 0, offset_b = 0;
						bool const cross = libmaus2::lcs::AlignmentTraceContainer::cross(
							ATC_a,start_a_r /* ref pos */,start_a_q /* read pos */,offset_a,
							ATC_b,start_b_r /* ref pos */,start_b_q /* read pos */,offset_b
						);

						int64_t const score_a = libmaus2::lcs::AlignmentTraceContainer::getScore(ATC_a.ta,ATC_a.te,1,1,1,1);
						int64_t const score_b = libmaus2::lcs::AlignmentTraceContainer::getScore(ATC_b.ta,ATC_b.te,1,1,1,1);

						if ( cross )
						{
							mindist = 0;

							if ( bestdist >= 0 )
							{
								best = it;
								bestdist = -1;
							}
						}
						else
						{
							libmaus2::math::IntegerInterval<int64_t> const aintv = prim_a->getReferenceInterval();
							libmaus2::math::IntegerInterval<int64_t> const bintv = prim_b->getReferenceInterval();

							if ( ! aintv.intersection(bintv).isEmpty() )
							{
								if ( bestdist > 0 )
								{
									mindist = 0;
									bestdist = 0;
									best = it;
								}
							}
							else
							{
								int64_t avdist;

								if ( aintv.from > bintv.to )
									avdist = aintv.from- bintv.to;
								else
									avdist = bintv.from - aintv.to;

								if ( avdist < mindist )
								{
									mindist = avdist;
									bestdist = avdist;
									best = it;
								}
							}
						}

						anycross = anycross || cross;
					}
					#endif
				}

				if ( anycross )
				{
					numcross += 1;

					uint64_t const a_size = Amappos.size();
					uint64_t const ab_size = ABmappos.size();

					double const frac = a_size ? static_cast<double>(ab_size)/static_cast<double>(a_size) : 0.0;
					fracsum += frac;
					std::cerr << "[V] frac " << frac << std::endl;
				}
				else
				{
					std::cerr << "[E] no cross " << (*a_itc)->getReferenceLength() << " correct seq " << numcorrectseq
						<< " erate " << (*a_itc)->getErrorRate() << " maps " << (b_top-b_itc) << " mindist " << mindist << " best dist " << bestdist
						<< (*a_itc)->getReferenceInterval()
						<< " "
						<< ((bestdist != std::numeric_limits<int64_t>::max()) ? (*best)->getReferenceInterval() : libmaus2::math::IntegerInterval<int64_t>::empty())
						<< " "
						<< (*a_itc)->getFrontSoftClipping() << " "
						<< ((bestdist != std::numeric_limits<int64_t>::max()) ? (*best)->getFrontSoftClipping() : 0)
						<< " "
						<< ((bestdist != std::numeric_limits<int64_t>::max()) ? (*best)->getCigarString() : std::string())
						<< std::endl;

					std::cerr << (*a_itc)->getCigarString() << std::endl;

					#if 0
					if ( (*a_itc)->getReferenceLength() < 2048 )
					{
						std::cerr << (*a_itc)->getRead() << std::endl;
					}
					#endif

					// if ( bestdist != std::numeric_limits<int64_t>::max() )
					if ( bestdist <= 0 )
					{
						std::vector<std::pair<int64_t,int64_t> > MPPAV, MPPBV;
						(*a_itc)->getMappingPositionPairs(MPPAV);
						(*best)->getMappingPositionPairs(MPPBV);

						uint64_t i = 0, j = 0;

						while ( i < MPPAV.size() && j < MPPBV.size() )
						{
							if ( MPPAV[i].first < MPPBV[j].first )
							{
								std::cerr << MPPAV[i] << std::endl;
								i++;
							}
							else if ( MPPAV[i].first > MPPBV[j].first )
							{
								std::cerr << "\t" << MPPBV[j] << std::endl;
								j++;
							}
							else
							{
								std::cerr << MPPAV[i] << "\t" << MPPBV[j] << std::endl;
								i++;
								j++;

							}
						}
						while ( i < MPPAV.size() )
						{
							std::cerr << MPPAV[i] << std::endl;
							i++;
						}
						while ( j < MPPBV.size() )
						{
							std::cerr << "\t" << MPPBV[j] << std::endl;
							j++;
						}
					}

					numnocross += 1;
				}

				#if 0
				libmaus2::bambam::BamAlignment::shared_ptr_type prim_a;
				libmaus2::bambam::BamAlignment::shared_ptr_type prim_b;

				for ( iterator it = a_itc; it != a_top; ++it )
					if ( ! (*it)->isSecondary() )
					{
						prim_a = *it;
						break;
					}

				for ( iterator it = b_itc; it != b_top; ++it )
					if ( ! (*it)->isSecondary() )
					{
						prim_b = *it;
						break;
					}

				if ( prim_a && prim_b )
				{
					#if 0
					prim_a->formatAlignment(std::cerr,header_a,aux);
					std::cerr.put('\n');
					prim_b->formatAlignment(std::cerr,header_b,aux);
					std::cerr.put('\n');
					#endif

					if ( prim_a->getRefID() != prim_b->getRefID() )
					{
						std::cerr << "[D] ref id mismatch" << std::endl;
						numwrongrefid += 1;
					}
					else if ( prim_a->isReverse() != prim_b->isReverse() )
					{
						std::cerr << "[D] strand mismatch" << std::endl;
						numwrongstrand += 1;
					}
					else
					{
						uint64_t const start_a_r = prim_a->getPos();
						uint64_t const start_a_q = keep_soft_clipping ? prim_a->getFrontSoftClipping() : 0;
						uint64_t const start_b_r = prim_b->getPos();
						uint64_t const start_b_q = keep_soft_clipping ? prim_b->getFrontSoftClipping() : 0;

						uint32_t const num_cig_a = prim_a->getCigarOperations(Acigop);
						uint32_t const num_cig_b = prim_b->getCigarOperations(Bcigop);

						libmaus2::bambam::CigarStringParser::cigarToTrace(Acigop.begin(),Acigop.begin()+num_cig_a,ATC_a,true /* ignore unknown */);
						libmaus2::bambam::CigarStringParser::cigarToTrace(Bcigop.begin(),Bcigop.begin()+num_cig_b,ATC_b,true /* ignore unknown */);

						uint64_t offset_a = 0, offset_b = 0;
						bool const cross = libmaus2::lcs::AlignmentTraceContainer::cross(
							ATC_a,start_a_r,start_a_q,offset_a,
							ATC_b,start_b_r,start_b_q,offset_b
						);

						int64_t const score_a = libmaus2::lcs::AlignmentTraceContainer::getScore(ATC_a.ta,ATC_a.te,1,1,1,1);
						int64_t const score_b = libmaus2::lcs::AlignmentTraceContainer::getScore(ATC_b.ta,ATC_b.te,1,1,1,1);

						if ( cross )
						{
							std::cerr << "[D] crossing " << prim_a->getReferenceLength() << " " << prim_b->getReferenceLength() << " " << score_a << " " << score_b << std::endl;
							numcross += 1;
						}
						else
						{
							libmaus2::math::IntegerInterval<int64_t> IA(prim_a->getPos(),prim_a->getPos()+prim_a->getReferenceLength()-1);
							libmaus2::math::IntegerInterval<int64_t> IB(prim_b->getPos(),prim_b->getPos()+prim_b->getReferenceLength()-1);
							std::cerr << "[D] no crossing " << IA << "," << IB << "," << prim_a->getFrontSoftClipping() << "," << prim_b->getFrontSoftClipping() << std::endl;
							numnocross += 1;
						}
					}
				}
				else if ( ! prim_a )
				{
					std::cerr << "no prim a" << std::endl;
					prim_b->formatAlignment(std::cerr,header_b,aux);
					std::cerr.put('\n');
					numnoprima += 1;
				}
				else if ( ! prim_b )
				{
					std::cerr << "no prim_b" << std::endl;
					prim_a->formatAlignment(std::cerr,header_a,aux);
					std::cerr.put('\n');
					numnoprimb += 1;
				}
				#endif

				a_itc = a_top;
				b_itc = b_top;
			}
			else
			{
				int const cmp = strcmpnum.strcmpnum((*a_itc)->getName(),(*b_itc)->getName());

				if ( cmp < 0 )
				{
					++a_itc;
					++skipa;
				}
				else if ( cmp > 0 )
				{
					++b_itc;
					++skipb;
				}
			}
		}

		std::cerr << "[V] numcross=" << numcross << std::endl;
		std::cerr << "[V] numnocross=" << numnocross << std::endl;
		std::cerr << "[V] numwrongrefid=" << numwrongrefid << std::endl;
		std::cerr << "[V] numwrongstrand=" << numwrongstrand << std::endl;
		std::cerr << "[V] skipa=" << skipa << std::endl;
		std::cerr << "[V] skipb=" << skipb << std::endl;
		std::cerr << "[V] avg frac=" << (numcross ? (fracsum / numcross) : 0.0) << std::endl;
		#endif
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
