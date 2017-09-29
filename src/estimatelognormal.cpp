/*
    loresim2
    Copyright (C) 2017 German Tischler

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
#include <libmaus2/random/LogNormalRandom.hpp>
#include <libmaus2/fastx/StreamFastAReader.hpp>
#include <iostream>

int main()
{
	try
	{
		std::vector<uint64_t> V;

		libmaus2::fastx::StreamFastAReaderWrapper SFAR(std::cin);
		libmaus2::fastx::StreamFastAReaderWrapper::pattern_type pattern;
		uint64_t tl = 0;

		while ( SFAR.getNextPatternUnlocked(pattern) )
		{
			uint64_t const l = pattern.spattern.size();

			if ( l )
			{
				while ( ! (l < V.size()) )
					V.push_back(0);

				assert ( l < V.size() );

				V[l]++;
				tl += l;
			}
		}

		std::pair<double,double> const P = libmaus2::random::LogNormalRandom::computeAverageAndSigma(V);
		std::pair<double,double> const PV = libmaus2::random::LogNormalRandom::computeParameters(P.first,P.second);
		std::cout << PV.first << "\t" << PV.second << "\t" << P.first << "\t" << P.second << std::endl;
	}
	catch(std::exception const & ex)
	{
		std::cerr << ex.what() << std::endl;
		return EXIT_FAILURE;
	}
}
