/*
    loresim
    Copyright (C) 2016 German Tischler

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
#if ! defined(INTV_HPP)
#define INTV_HPP

#include <libmaus2/math/IntegerInterval.hpp>
#include <libmaus2/util/stringFunctions.hpp>
#include <libmaus2/bambam/BamAlignment.hpp>
#include <sstream>

struct Intv
{
	bool valid;
	libmaus2::math::IntegerInterval<int64_t> intv;
	double erate;
	int64_t numerr;
	int64_t ins;
	int64_t sub;
	int64_t del;

	static bool expect(std::istream & in, std::string const & s)
	{
		uint64_t i = 0;
		bool ok = true;
		while ( ok && i < s.size() && in.peek() != std::istream::traits_type::eof() )
			ok = ok && (s[i++] == in.get());
		return ok;
	}

	static bool getNumber(std::istream & in, int64_t & num)
	{
		bool ok = false;
		num = 0;
		while ( in.peek() != std::istream::traits_type::eof() && isdigit(in.peek()) )
		{
			num *= 10;
			num += in.get()-'0';
			ok = true;
		}

		return ok;
	}

	static bool getUntilTerm(std::istream & in, std::ostream & ostr, int term)
	{
		while ( in.peek() != std::istream::traits_type::eof() && in.peek() != term )
			ostr.put(in.get());

		if ( in.peek() != std::istream::traits_type::eof() )
		{
			int const t = in.get();
			assert ( t == term );
			return true;
		}
		else
		{
			return false;
		}
	}

	static bool getFloat(std::istream & in, double & d, int term)
	{
		std::ostringstream ostr;
		bool ok = getUntilTerm(in,ostr,term);
		std::istringstream istr(ostr.str());
		istr >> d;
		ok = ok && istr;
		return ok;
	}

	bool setup(std::istream & in)
	{
		bool ok = true;

		ok = ok && expect(in,"intv([");
		ok = ok && getNumber(in,intv.from);
		ok = ok && expect(in,",");
		ok = ok && getNumber(in,intv.to);
		if ( ok )
			intv.to -= 1;
		ok = ok && expect(in,"),erate=");
		ok = ok && getFloat(in,erate,',');
		ok = ok && expect(in,"numerr=");
		ok = ok && getNumber(in,numerr);
		ok = ok && expect(in,",ins=");
		ok = ok && getNumber(in,ins);
		ok = ok && expect(in,",sub=");
		ok = ok && getNumber(in,sub);
		ok = ok && expect(in,",del=");
		ok = ok && getNumber(in,del);
		ok = ok && expect(in,")");
		return ok;
	}

	Intv() : valid(false)
	{

	}

	Intv(std::istream & in)
	: valid(false)
	{
		valid = setup(in);
	}

	Intv(std::string const & s)
	: valid(false)
	{
		std::istringstream istr(s);
		valid = setup(istr);
	}

	static std::vector<Intv> parse(std::string const & s)
	{
		std::vector<Intv> Vintv;

		std::deque<std::string> tokens = libmaus2::util::stringFunctions::tokenize<std::string>(s,std::string(";"));

		for ( uint64_t i = 0; i < tokens.size(); ++i )
			if ( tokens[i].size() && tokens[i].size() >= strlen("intv") && tokens[i].substr(0,strlen("intv")) == std::string("intv") )
			{
				Intv intv(tokens[i]);
				if (intv.valid )
					Vintv.push_back(intv);
			}

		for ( uint64_t i = 1; i < Vintv.size(); ++i )
			assert ( Vintv[i-1].intv.to+1 == Vintv[i].intv.from );

		return Vintv;
	}

	static std::vector<Intv> computeRIntv(std::vector<Intv> const & Vintv, libmaus2::bambam::BamAlignment const & a_algn)
	{
		return computeRIntv(Vintv, a_algn.D.get());
	}

	static std::vector<Intv> computeRIntv(std::vector<Intv> const & Vintv, uint8_t const * D)
	{
		assert ( static_cast<int64_t>(libmaus2::bambam::BamAlignmentDecoderBase::getReferenceLength(D)) == Vintv.back().intv.to+1 );

		libmaus2::autoarray::AutoArray<libmaus2::bambam::cigar_operation> cigop;
		uint64_t const numcig = libmaus2::bambam::BamAlignmentDecoderBase::getCigarOperations(D,cigop);
		std::vector < libmaus2::bambam::BamFlagBase::bam_cigar_ops > expcig;
		for ( uint64_t i = 0; i < numcig; ++i )
			for ( int64_t j = 0; j < cigop[i].second; ++j )
				expcig.push_back(static_cast<libmaus2::bambam::BamFlagBase::bam_cigar_ops>(cigop[i].first));
		uint64_t zz = 0;
		while (
			zz < expcig.size() &&
			(expcig[zz] == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CSOFT_CLIP ||
			 expcig[zz] == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CHARD_CLIP )
		)
			++zz;

		std::vector<Intv> Rintv;
		uint64_t rlsum = 0;
		for ( uint64_t zi = 0; zi < Vintv.size(); ++zi )
		{
			// length on reference
			uint64_t d = Vintv[zi].intv.to-Vintv[zi].intv.from+1;
			// length on read
			uint64_t rl = 0;

			uint64_t ins = 0, del = 0, match = 0, mismatch = 0;

			while ( d )
			{
				switch ( expcig[zz++] )
				{
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CMATCH:
						--d;
						++rl;
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CEQUAL:
						match += 1;
						--d;
						++rl;
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDIFF:
						mismatch += 1;
						--d;
						++rl;
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS:
						ins += 1;
						++rl;
						break;
					case libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CDEL:
						del += 1;
						--d;
						break;
					default:
						break;
				}
			}

			while ( zz < expcig.size() && expcig[zz] == libmaus2::bambam::BamFlagBase::LIBMAUS2_BAMBAM_CINS )
			{
				++zz;
				++rl;
			}

			Intv rintv;
			rintv.valid = true;
			rintv.intv.from = rlsum;
			rintv.intv.to = rintv.intv.from + rl - 1;
			rintv.ins = ins;
			rintv.del = del;
			rintv.sub = mismatch;
			rintv.numerr = ins+del+mismatch;
			rintv.erate = (ins+del+mismatch) / static_cast<double>(ins+del+mismatch+match);
			Rintv.push_back(rintv);

			rlsum += rl;
		}

		assert (
			static_cast<int64_t>(
				rlsum+
				libmaus2::bambam::BamAlignmentDecoderBase::getFrontSoftClipping(D) +
				libmaus2::bambam::BamAlignmentDecoderBase::getBackSoftClipping(D)
			)
			==
			libmaus2::bambam::BamAlignmentDecoderBase::getLseq(D)
		);

		return Rintv;
	}
};

inline std::ostream & operator<<(std::ostream & out, Intv const & intv)
{
	return out << "intv([" << intv.intv.from << "," << intv.intv.to+1 << "),erate=" << intv.erate << ",numerr=" << intv.numerr << ",ins=" << intv.ins << ",sub=" << intv.sub << ",del=" << intv.del << ")";
}
#endif
