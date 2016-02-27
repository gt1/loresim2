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
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <libmaus2/fastx/StreamFastQReader.hpp>
#include <libmaus2/aio/PosixFdInputStream.hpp>
#include <libmaus2/fastx/SpaceTable.hpp>
#include <libmaus2/random/Random.hpp>

template<typename reader_type>
void fastareformat(libmaus2::util::ArgInfo const & arginfo, reader_type & reader)
{
	typename reader_type::pattern_type pattern;

	std::string const prolog = arginfo.getUnparsedValue("prolog","L");
	uint64_t readid = arginfo.getValueUnsignedNumeric<uint64_t>("readidbase",0);
	uint64_t const cols = arginfo.getValueUnsignedNumeric<uint64_t>("cols",80);
	bool const randomN = arginfo.getValueUnsignedNumeric<uint64_t>("randomN",true);
	libmaus2::fastx::SpaceTable const ST;
	libmaus2::random::Random::setup();
	bool replaced = false;
	bool replaceprinted = false;

	while ( reader.getNextPatternUnlocked(pattern) )
	{
		std::string s = pattern.spattern;

		if ( randomN )
			for ( uint64_t i = 0; i < s.size(); ++i )
				if ( libmaus2::fastx::mapChar(s[i]) >= 4 )
				{
					s[i] = libmaus2::fastx::remapChar(libmaus2::random::Random::rand8() & 3);
					replaced = true;
				}
		char const * c = s.c_str();

		if ( replaced && ! replaceprinted )
		{
			std::cerr << "[V] warning, replacing non ACGT symbols with random bases" << std::endl;
			replaceprinted = true;
		}

		std::cout << '>' << prolog << '/' << (readid++) << '/' << 0 << '_' << s.size() << " RQ=0.851 " << pattern.sid << "\n";
		uint64_t low = 0;

		while ( low < s.size() )
		{
			uint64_t const high = std::min(low+cols,static_cast<uint64_t>(s.size()));

			std::cout.write(c+low,high-low);
			std::cout.put('\n');

			low = high;
		}
	}
}

int main(int argc, char * argv[])
{
	try
	{
		libmaus2::util::ArgInfo const arginfo(argc,argv);

		libmaus2::aio::PosixFdInputStream PFIS(STDIN_FILENO,64*1024,1024);

		int const c = PFIS.peek();

		if ( c == '>' || c == std::istream::traits_type::eof() )
		{
			libmaus2::fastx::StreamFastAReaderWrapper in(PFIS);
			fastareformat(arginfo,in);
		}
		else if ( c == '@' )
		{
			libmaus2::fastx::StreamFastQReaderWrapper in(PFIS);
			fastareformat(arginfo,in);
		}
		else
		{
			libmaus2::exception::LibMausException lme;
			lme.getStream() << "Unknown input file format, first character " << static_cast<char>(c) << std::endl;
			lme.finish();
			throw lme;
		}
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
