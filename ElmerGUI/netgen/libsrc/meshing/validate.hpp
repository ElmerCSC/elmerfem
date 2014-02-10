#ifndef VALIDATE_HPP
#define VALIDATE_HPP



void GetPureBadness(Mesh & mesh, ARRAY<double> & pure_badness,
		    const BitArray & isnewpoint);
double Validate(const Mesh & mesh, ARRAY<ElementIndex> & bad_elements,
		const ARRAY<double> & pure_badness, 
		double max_worsening, const bool uselocalworsening,
		ARRAY<double> * quality_loss = NULL);
void RepairBisection(Mesh & mesh, ARRAY<ElementIndex> & bad_elements, const BitArray & isnewpoint, Refinement & refinement,
		     const ARRAY<double> & pure_badness, 
		     double max_worsening, const bool uselocalworsening,
		     const ARRAY< ARRAY<int,PointIndex::BASE>* > & idmaps);

#endif // VALIDATE_HPP
