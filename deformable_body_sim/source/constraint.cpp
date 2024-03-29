// ---------------------------------------------------------------------------------//
// Copyright (c) 2013, Regents of the University of Pennsylvania                    //
// All rights reserved.                                                             //
//                                                                                  //
// Redistribution and use in source and binary forms, with or without               //
// modification, are permitted provided that the following conditions are met:      //
//     * Redistributions of source code must retain the above copyright             //
//       notice, this list of conditions and the following disclaimer.              //
//     * Redistributions in binary form must reproduce the above copyright          //
//       notice, this list of conditions and the following disclaimer in the        //
//       documentation and/or other materials provided with the distribution.       //
//     * Neither the name of the <organization> nor the                             //
//       names of its contributors may be used to endorse or promote products       //
//       derived from this software without specific prior written permission.      //
//                                                                                  //
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  //
// ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    //
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           //
// DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY               //
// DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       //
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     //
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      //
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       //
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    //
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     //
//                                                                                  //
// Contact Tiantian Liu (ltt1598@gmail.com) if you have any questions.              //
//----------------------------------------------------------------------------------//

#include "constraint.h"

//----------Constraint Class----------//
Constraint::Constraint(ScalarType *stiffness, ScalarType *pbd_stiffness) : 
    m_p_stiffness(stiffness),
	m_p_pbd_stiffness(pbd_stiffness)
{
}

Constraint::Constraint(const Constraint& other) : 
    m_p_stiffness(other.m_p_stiffness),
	m_p_pbd_stiffness(other.m_p_pbd_stiffness)
{
}

Constraint::~Constraint()
{
}

//----------AttachmentConstraint Class----------//
AttachmentConstraint::AttachmentConstraint(ScalarType *stiffness, ScalarType *pbd_stiffness) : 
    Constraint(stiffness, pbd_stiffness)
{
    m_selected = false;
}

AttachmentConstraint::AttachmentConstraint(ScalarType *stiffness, ScalarType *pbd_stiffness, unsigned int p0, const EigenVector3& fixedpoint) : 
    Constraint(stiffness, pbd_stiffness),
    m_p0(p0),
    m_fixd_point(fixedpoint)
{
    m_selected = false;
}

AttachmentConstraint::AttachmentConstraint(const AttachmentConstraint& other) : 
    Constraint(other),
    m_p0(other.m_p0),
    m_fixd_point(other.m_fixd_point),
    m_selected(other.m_selected)
{
    
}

AttachmentConstraint::~AttachmentConstraint()
{
}

// attachment spring gradient: k*(current_length)*current_direction
void AttachmentConstraint::EvaluateGradient(const VectorX& x, VectorX& gradient)
{
	// LOOK
    EigenVector3 g_i = (*(m_p_stiffness))*(x.block_vector(m_p0) - m_fixd_point);
    gradient.block_vector(m_p0) += g_i;
}

void AttachmentConstraint::EvaluateHessian(const VectorX& x, std::vector<SparseMatrixTriplet>& hessian_triplets)
{
	// LOOK
    ScalarType ks = *(m_p_stiffness);
    hessian_triplets.push_back(SparseMatrixTriplet(3*m_p0, 3*m_p0, ks));
    hessian_triplets.push_back(SparseMatrixTriplet(3*m_p0+1, 3*m_p0+1, ks));
    hessian_triplets.push_back(SparseMatrixTriplet(3*m_p0+2, 3*m_p0+2, ks));
}

void AttachmentConstraint::PBDProject(VectorX& x, const SparseMatrix& inv_mass, unsigned int ns)
{
	// LOOK
	ScalarType k_prime = 1 - std::pow(1-*(m_p_pbd_stiffness), 1.0/ns);
	
	EigenVector3 p = x.block_vector(m_p0);
	EigenVector3 dp = m_fixd_point-p;

	x.block_vector(m_p0) += k_prime * dp;
}

void AttachmentConstraint::Draw(const VBO& vbos)
{
    m_attachment_constraint_body.move_to(Eigen2GLM(m_fixd_point));
    if (m_selected)
        m_attachment_constraint_body.change_color(glm::vec3(0.8, 0.8, 0.2));
    else
        m_attachment_constraint_body.change_color(glm::vec3(0.8, 0.2, 0.2));
        
    m_attachment_constraint_body.Draw(vbos);
}

//----------SpringConstraint Class----------//
SpringConstraint::SpringConstraint(ScalarType *stiffness, ScalarType *pbd_stiffness) : 
    Constraint(stiffness, pbd_stiffness)
{
}

SpringConstraint::SpringConstraint(ScalarType *stiffness, ScalarType *pbd_stiffness, unsigned int p1, unsigned int p2, ScalarType length) : 
    Constraint(stiffness, pbd_stiffness),
    m_p1(p1),
    m_p2(p2),
    m_rest_length(length)
{
}

SpringConstraint::SpringConstraint(const SpringConstraint& other) : 
    Constraint(other),
    m_p1(other.m_p1),
    m_p2(other.m_p2),
    m_rest_length(other.m_rest_length)
{
}

SpringConstraint::~SpringConstraint()
{
}

// sping gradient: k*(current_length-rest_length)*current_direction;
void SpringConstraint::EvaluateGradient(const VectorX& x, VectorX& gradient)
{
	// TODO
	EigenVector3 vec_x_ij = x.block_vector(m_p1) - x.block_vector(m_p2);
	
	double k = *(m_p_stiffness);
	double current_length = vec_x_ij.norm();
	EigenVector3 current_direction = vec_x_ij.normalized();
	EigenVector3 gradient_ij = k * (current_length - m_rest_length) * current_direction;
	gradient.block_vector(m_p1) += gradient_ij;
	gradient.block_vector(m_p2) -= gradient_ij;
}

void SpringConstraint::EvaluateHessian(const VectorX& x, std::vector<SparseMatrixTriplet>& hessian_triplets)
{
	// TODO
	EigenVector3 x_ij = x.block_vector(m_p1) - x.block_vector(m_p2);
    double l_ij = x_ij.norm();
	double ratio =  m_rest_length/l_ij;

    EigenMatrix3 k = *(m_p_stiffness)*(EigenMatrix3::Identity() - ratio*(EigenMatrix3::Identity() - (x_ij*x_ij.transpose())/pow(l_ij, 2)));

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            ScalarType val = k(i,j);
            hessian_triplets.push_back(SparseMatrixTriplet(3*m_p1+i, 3*m_p1+j, val));
			hessian_triplets.push_back(SparseMatrixTriplet(3*m_p2+i, 3*m_p2+j, val));
            hessian_triplets.push_back(SparseMatrixTriplet(3*m_p1+i, 3*m_p2+j, -val));
            hessian_triplets.push_back(SparseMatrixTriplet(3*m_p2+i, 3*m_p1+j, -val));
        }
    }

}

void SpringConstraint::PBDProject(VectorX& x, const SparseMatrix& inv_mass, unsigned int ns)
{
	// TODO

	ScalarType k_prime = 1 - std::pow(1-*(m_p_pbd_stiffness), 1.0/ns);
	
	EigenVector3 p1 = x.block_vector(m_p1);
	EigenVector3 p2 = x.block_vector(m_p2);

	EigenVector3 dp = p1-p2;
	double dp_dist = dp.norm();

	double w1 = inv_mass.coeff(m_p1, m_p1);
	double w2 = inv_mass.coeff(m_p2, m_p2);

	if(dp.norm() < 1e-3){
		EigenVector3 delta_p1 = (-w1/(w1+w2)) * (dp_dist + 1e-3 - m_rest_length) * dp/(dp_dist + 1e-3);
		EigenVector3 delta_p2 = (w1/(w1+w2)) * (dp_dist  + 1e-3 - m_rest_length) * dp/(dp_dist + 1e-3);

		x.block_vector(m_p1) += k_prime * delta_p1;
		x.block_vector(m_p2) += k_prime * delta_p2;
	}else{
		// divide by 0
		EigenVector3 delta_p1 = (-w1/(w1+w2)) * (dp_dist - m_rest_length) * dp.normalized();
		EigenVector3 delta_p2 = (w1/(w1+w2)) * (dp_dist - m_rest_length) * dp.normalized();

		x.block_vector(m_p1) += k_prime * delta_p1;
		x.block_vector(m_p2) += k_prime * delta_p2;
	}
}