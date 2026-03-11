/*
* Copyright (C) 2020-2026 MEmilio
*
* Authors: Daniel Abele, Jan Kleinert, Martin J. Kuehn
*
* Contact: Martin J. Kuehn <Martin.Kuehn@DLR.de>
*
* Licensed under the Apache License, Version 2.0 (the "License");
* you may not use this file except in compliance with the License.
* You may obtain a copy of the License at
*
*     http://www.apache.org/licenses/LICENSE-2.0
*
* Unless required by applicable law or agreed to in writing, software
* distributed under the License is distributed on an "AS IS" BASIS,
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
* See the License for the specific language governing permissions and
* limitations under the License.
*/
#ifndef SEIR_VECTOR_PARAMETERS_H
#define SEIR_VECTOR_PARAMETERS_H

#include "memilio/config.h"
#include "memilio/epidemiology/age_group.h"
#include "memilio/epidemiology/uncertain_matrix.h"
#include "memilio/utils/custom_index_array.h"
#include "memilio/utils/uncertain_value.h"
#include "memilio/utils/parameter_set.h"

namespace mio
{
namespace oseirvector
{

/***************************************
 * Define Parameters of the SEIR model *
 ***************************************/

/**
 * @brief probability of getting infected from a contact
 */
template <typename FP>
struct TransmissionProbabilityOnContact {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 1.);
    }
    static std::string name()
    {
        return "TransmissionProbabilityOnContact";
    }
};

/**
 * @brief the latent time in day unit
 */
template <typename FP>
struct TimeExposed {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 5.2);
    }
    static std::string name()
    {
        return "TimeExposed";
    }
};

/**
 * @brief the infectious time in day unit
 */
template <typename FP>
struct TimeInfected {
    using Type = CustomIndexArray<UncertainValue<FP>, AgeGroup>;
    static Type get_default(AgeGroup size)
    {
        return Type(size, 6.0);
    }
    static std::string name()
    {
        return "TimeInfected";
    }
};

/**
 * @brief the contact patterns within the society are modelled using a ContactMatrix
 */
template <class FP>
struct ContactPatterns {
    using Type = UncertainContactMatrix<FP>;
    static Type get_default(AgeGroup size)
    {
        return Type(1, static_cast<Eigen::Index>((size_t)size));
    }
    static std::string name()
    {
        return "ContactPatterns";
    }
};

// New Parameters added 

/**
 * @brief birth rate of the mosquito
 */
template <typename FP>
struct MosquitoBirthRate {
    using Type = UncertainValue<FP>;
    static Type get_default() { return Type(0.1); }  // example default value
    static std::string name() { return "MosquitoBirthRate"; }
};


/**
 * @brief death rate of the mosquito
 */
template <typename FP>
struct MosquitoDeathRate {
    using Type = UncertainValue<FP>;
    static Type get_default() { return Type(0.1); }  // example default value assuming btw that it is equal to the birth rate
    static std::string name() { return "MosquitoDeathRate"; }
};

/**
 * @brief Transmission from Humans to Vector
 */
template <typename FP>
struct TransmissionHumanToVector {
    using Type = UncertainValue<FP>;
    static Type get_default() { return Type(0.3); }  // example default value 
    static std::string name() { return "TransmissionHumanToVector"; }
};

/**
 * @brief Transmission from Vectors to Humans
 */
template <typename FP>
struct TransmissionVectorToHuman {
    using Type = UncertainValue<FP>;
    static Type get_default() { return Type(0.2); }  // example default value 
    static std::string name() { return "TransmissionVectorToHuman"; }
};

/**
 * @brief Mosquito biting rate
 */
template <typename FP>
struct MosquitoBitingRate {
    using Type = UncertainValue<FP>;
    static Type get_default() { return Type(0.5); }  // example default value 
    static std::string name() { return "MosquitoBitingRate"; }
};


template <typename FP>
using ParametersBase =
    ParameterSet<TransmissionProbabilityOnContact<FP>, TimeExposed<FP>, TimeInfected<FP>, ContactPatterns<FP>,
     MosquitoBirthRate<FP>, TransmissionHumanToVector<FP>, TransmissionVectorToHuman<FP>, MosquitoDeathRate<FP>, 
     MosquitoBitingRate<FP>>;

/**
 * @brief Parameters of an age-resolved SECIR/SECIHURD model.
 */
template <typename FP>
class Parameters : public ParametersBase<FP>
{
public:
    Parameters(AgeGroup num_agegroups)
        : ParametersBase<FP>(num_agegroups)
        , m_num_groups{num_agegroups}
    {
    }

    AgeGroup get_num_groups() const
    {
        return m_num_groups;
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and applies them, if they do not.
     * Time spans cannot be negative and probabilities can only take values between [0,1].
     *
     * Attention: This function should be used with care. It is necessary for some test problems to run through quickly,
     *            but in a manual execution of an example, check_constraints() may be preferred. Note that the apply_constraints()
     *            function can and will not set Parameters to meaningful values in an epidemiological or virological context,
     *            as all models are designed to be transferable to multiple diseases. Consequently, only acceptable
     *            (like 0 or 1 for probabilities or small positive values for time spans) values are set here and a manual adaptation
     *            may often be necessary to have set meaningful values.
     *
     * @return Returns true if one ore more constraint were corrected, false otherwise.
     */
    bool apply_constraints()
    {
        const FP tol_times = 1e-1;
        const FP tol_pos   = 1e-6;
        bool corrected = false;

        for (auto i = AgeGroup(0); i < AgeGroup(m_num_groups); ++i) {
            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeExposed changed from {} to {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeExposed<FP>>()[i], tol_times);
                this->template get<TimeExposed<FP>>()[i] = tol_times;
                corrected                                = true;
            }
            if (this->template get<TimeInfected<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfected changed from {} to {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfected<FP>>()[i], tol_times);
                this->template get<TimeInfected<FP>>()[i] = tol_times;
                corrected                                 = true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_warning("Constraint check: Parameter TransmissionProbabilityOnContact changed from {} to {} ",
                            this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0.0);
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] = 0.0;
                corrected                                                     = true;
            }
        }
        // New paramters mosquito related 
            if (this->template get<MosquitoBirthRate<FP>>() < 0.0) {
            log_warning("Constraint check: MosquitoBirthRate {} smaller than 0. Setting to 0.",
                        this->template get<MosquitoBirthRate<FP>>());
            this->template get<MosquitoBirthRate<FP>>() = 0.0;
            corrected = true;
        }

        if (this->template get<MosquitoDeathRate<FP>>() <= 0.0) {
            log_warning("Constraint check: MosquitoDeathRate {} must be > 0. Setting to {}.",
                        this->template get<MosquitoDeathRate<FP>>(), tol_pos);
            this->template get<MosquitoDeathRate<FP>>() = tol_pos;
            corrected = true;
        }

        if (this->template get<MosquitoBitingRate<FP>>() <= 0.0) {
            log_warning("Constraint check: MosquitoBitingRate {} must be > 0. Setting to {}.",
                        this->template get<MosquitoBitingRate<FP>>(), tol_pos);
            this->template get<MosquitoBitingRate<FP>>() = tol_pos;
            corrected = true;
        }

        if (this->template get<TransmissionVectorToHuman<FP>>() < 0.0 ||
            this->template get<TransmissionVectorToHuman<FP>>() > 1.0) {
            log_warning("Constraint check: TransmissionVectorToHuman {} outside [0,1]. Setting to 0.",
                        this->template get<TransmissionVectorToHuman<FP>>());
            this->template get<TransmissionVectorToHuman<FP>>() = 0.0;
            corrected = true;
        }

        if (this->template get<TransmissionHumanToVector<FP>>() < 0.0 ||
            this->template get<TransmissionHumanToVector<FP>>() > 1.0) {
            log_warning("Constraint check: TransmissionHumanToVector {} outside [0,1]. Setting to 0.",
                        this->template get<TransmissionHumanToVector<FP>>());
            this->template get<TransmissionHumanToVector<FP>>() = 0.0;
            corrected = true;
        }






        return corrected;
    }

    /**
     * @brief Checks whether all Parameters satisfy their corresponding constraints and logs an error
     * if constraints are not satisfied.
     * @return Returns true if one constraint is not satisfied, otherwise false.
     */
    bool check_constraints() const
    {
        const FP tol_times = 1e-1;
       // const FP tol_pos   = 1e-6;

        for (auto i = AgeGroup(0); i < m_num_groups; i++) {
            if (this->template get<TimeExposed<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeExposed {} smaller or equal {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeExposed<FP>>()[i], tol_times);
                return true;
            }
            if (this->template get<TimeInfected<FP>>()[i] < tol_times) {
                log_warning(
                    "Constraint check: Parameter TimeInfected {} smaller or equal {}. Please note that "
                    "unreasonably small compartment stays lead to massively increased run time. Consider to cancel "
                    "and reset parameters.",
                    this->template get<TimeInfected<FP>>()[i], tol_times);
                return true;
            }
            if (this->template get<TransmissionProbabilityOnContact<FP>>()[i] < 0.0 ||
                this->template get<TransmissionProbabilityOnContact<FP>>()[i] > 1.0) {
                log_error("Constraint check: Parameter TransmissionProbabilityOnContact {} smaller {} or "
                          "greater {}",
                          this->template get<TransmissionProbabilityOnContact<FP>>()[i], 0.0, 1.0);
                return true;
            }
        }

        // New parameters
            if (this->template get<MosquitoBirthRate<FP>>() < 0.0) {
            log_error("Constraint check: MosquitoBirthRate {} must be >= 0.",
                    this->template get<MosquitoBirthRate<FP>>());
            return true;
        }

        if (this->template get<MosquitoDeathRate<FP>>() <= 0.0) {
            log_error("Constraint check: MosquitoDeathRate {} must be > 0.",
                    this->template get<MosquitoDeathRate<FP>>());
            return true;
        }

        if (this->template get<MosquitoBitingRate<FP>>() <= 0.0) {
            log_error("Constraint check: MosquitoBitingRate {} must be > 0.",
                    this->template get<MosquitoBitingRate<FP>>());
            return true;
        }

        if (this->template get<TransmissionVectorToHuman<FP>>() < 0.0 ||
            this->template get<TransmissionVectorToHuman<FP>>() > 1.0) {
            log_error("Constraint check: TransmissionVectorToHuman {} outside [0,1].",
                    this->template get<TransmissionVectorToHuman<FP>>());
            return true;
        }

        if (this->template get<TransmissionHumanToVector<FP>>() < 0.0 ||
            this->template get<TransmissionHumanToVector<FP>>() > 1.0) {
            log_error("Constraint check: TransmissionHumanToVector {} outside [0,1].",
                    this->template get<TransmissionHumanToVector<FP>>());
            return true;
        }


        return false;
    }

private:
    Parameters(ParametersBase<FP>&& base)
        : ParametersBase<FP>(std::move(base))
        , m_num_groups(this->template get<ContactPatterns<FP>>().get_cont_freq_mat().get_num_groups())
    {
    }

public:
    /**
     * deserialize an object of this class.
     * @see mio::deserialize
     */
    template <class IOContext>
    static IOResult<Parameters> deserialize(IOContext& io)
    {
        BOOST_OUTCOME_TRY(auto&& base, ParametersBase<FP>::deserialize(io));
        return success(Parameters(std::move(base)));
    }

private:
    AgeGroup m_num_groups;
};
} // namespace oseirvector
} // namespace mio

#endif // SEIR_VECTOR_PARAMETERS_H
