// Generated by rstantools.  Do not edit by hand.

/*
    tsnet is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    tsnet is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with tsnet.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef MODELS_HPP
#define MODELS_HPP
#define STAN__SERVICES__COMMAND_HPP
#ifndef USE_STANC3
#define USE_STANC3
#endif
#include <rstan/rstaninc.hpp>
// Code generated by stanc v2.32.2
#include <stan/model/model_header.hpp>
namespace model_VAR_wishart_namespace {
using stan::model::model_base_crtp;
using namespace stan::math;
stan::math::profile_map profiles__;
static constexpr std::array<const char*, 56> locations_array__ =
  {" (found before start of program)",
  " (in 'string', line 23, column 2 to column 23)",
  " (in 'string', line 28, column 2 to column 22)",
  " (in 'string', line 33, column 2 to column 67)",
  " (in 'string', line 35, column 2 to column 41)",
  " (in 'string', line 38, column 2 to column 18)",
  " (in 'string', line 71, column 2 to column 22)",
  " (in 'string', line 45, column 10 to column 23)",
  " (in 'string', line 44, column 13 to line 46, column 9)",
  " (in 'string', line 43, column 10 to column 65)",
  " (in 'string', line 42, column 18 to line 44, column 9)",
  " (in 'string', line 42, column 8 to line 46, column 9)",
  " (in 'string', line 41, column 19 to line 47, column 7)",
  " (in 'string', line 41, column 6 to line 47, column 7)",
  " (in 'string', line 40, column 17 to line 48, column 5)",
  " (in 'string', line 40, column 4 to line 48, column 5)",
  " (in 'string', line 39, column 2 to line 49, column 3)",
  " (in 'string', line 75, column 13 to column 14)",
  " (in 'string', line 75, column 6 to column 36)",
  " (in 'string', line 76, column 6 to column 59)",
  " (in 'string', line 73, column 17 to line 77, column 5)",
  " (in 'string', line 73, column 4 to line 77, column 5)",
  " (in 'string', line 72, column 2 to line 78, column 3)",
  " (in 'string', line 54, column 2 to column 50)",
  " (in 'string', line 60, column 2 to column 78)",
  " (in 'string', line 64, column 13 to column 14)",
  " (in 'string', line 64, column 6 to column 36)",
  " (in 'string', line 65, column 7 to column 54)",
  " (in 'string', line 62, column 17 to line 66, column 5)",
  " (in 'string', line 62, column 4 to line 66, column 5)",
  " (in 'string', line 61, column 2 to line 67, column 3)",
  " (in 'string', line 5, column 2 to column 17)",
  " (in 'string', line 6, column 2 to column 17)",
  " (in 'string', line 8, column 8 to column 9)",
  " (in 'string', line 8, column 18 to column 19)",
  " (in 'string', line 8, column 2 to column 23)",
  " (in 'string', line 10, column 9 to column 10)",
  " (in 'string', line 10, column 11 to column 12)",
  " (in 'string', line 10, column 2 to column 29)",
  " (in 'string', line 11, column 9 to column 10)",
  " (in 'string', line 11, column 11 to column 12)",
  " (in 'string', line 11, column 2 to column 31)",
  " (in 'string', line 14, column 2 to column 34)",
  " (in 'string', line 18, column 9 to column 10)",
  " (in 'string', line 18, column 11 to column 12)",
  " (in 'string', line 18, column 2 to column 48)",
  " (in 'string', line 23, column 9 to column 10)",
  " (in 'string', line 23, column 11 to column 12)",
  " (in 'string', line 28, column 13 to column 14)",
  " (in 'string', line 33, column 9 to column 10)",
  " (in 'string', line 33, column 11 to column 12)",
  " (in 'string', line 35, column 9 to column 10)",
  " (in 'string', line 35, column 11 to column 12)",
  " (in 'string', line 38, column 9 to column 10)",
  " (in 'string', line 38, column 11 to column 12)",
  " (in 'string', line 71, column 9 to column 12)"};
#include <stan_meta_header.hpp>
class model_VAR_wishart final : public model_base_crtp<model_VAR_wishart> {
private:
  int K;
  int T;
  std::vector<Eigen::Matrix<double,-1,1>> Y;
  Eigen::Matrix<double,-1,-1> prior_Beta_loc_data__;
  Eigen::Matrix<double,-1,-1> prior_Beta_scale_data__;
  int prior_Rho_marginal;
  Eigen::Matrix<double,-1,-1> I_data__;
  int log_lik_1dim__;
  Eigen::Map<Eigen::Matrix<double,-1,-1>> prior_Beta_loc{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> prior_Beta_scale{nullptr, 0, 0};
  Eigen::Map<Eigen::Matrix<double,-1,-1>> I{nullptr, 0, 0};
public:
  ~model_VAR_wishart() {}
  model_VAR_wishart(stan::io::var_context& context__, unsigned int
                    random_seed__ = 0, std::ostream* pstream__ = nullptr)
      : model_base_crtp(0) {
    int current_statement__ = 0;
    using local_scalar_t__ = double;
    boost::ecuyer1988 base_rng__ =
      stan::services::util::create_rng(random_seed__, 0);
    // suppress unused var warning
    (void) base_rng__;
    static constexpr const char* function__ =
      "model_VAR_wishart_namespace::model_VAR_wishart";
    // suppress unused var warning
    (void) function__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      current_statement__ = 31;
      context__.validate_dims("data initialization", "K", "int",
        std::vector<size_t>{});
      K = std::numeric_limits<int>::min();
      current_statement__ = 31;
      K = context__.vals_i("K")[(1 - 1)];
      current_statement__ = 31;
      stan::math::check_greater_or_equal(function__, "K", K, 0);
      current_statement__ = 32;
      context__.validate_dims("data initialization", "T", "int",
        std::vector<size_t>{});
      T = std::numeric_limits<int>::min();
      current_statement__ = 32;
      T = context__.vals_i("T")[(1 - 1)];
      current_statement__ = 32;
      stan::math::check_greater_or_equal(function__, "T", T, 0);
      current_statement__ = 33;
      stan::math::validate_non_negative_index("Y", "T", T);
      current_statement__ = 34;
      stan::math::validate_non_negative_index("Y", "K", K);
      current_statement__ = 35;
      context__.validate_dims("data initialization", "Y", "double",
        std::vector<size_t>{static_cast<size_t>(T), static_cast<size_t>(K)});
      Y = std::vector<Eigen::Matrix<double,-1,1>>(T,
            Eigen::Matrix<double,-1,1>::Constant(K,
              std::numeric_limits<double>::quiet_NaN()));
      {
        std::vector<local_scalar_t__> Y_flat__;
        current_statement__ = 35;
        Y_flat__ = context__.vals_r("Y");
        current_statement__ = 35;
        pos__ = 1;
        current_statement__ = 35;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 35;
          for (int sym2__ = 1; sym2__ <= T; ++sym2__) {
            current_statement__ = 35;
            stan::model::assign(Y, Y_flat__[(pos__ - 1)],
              "assigning variable Y", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 35;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 36;
      stan::math::validate_non_negative_index("prior_Beta_loc", "K", K);
      current_statement__ = 37;
      stan::math::validate_non_negative_index("prior_Beta_loc", "K", K);
      current_statement__ = 38;
      context__.validate_dims("data initialization", "prior_Beta_loc",
        "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(K)});
      prior_Beta_loc_data__ = Eigen::Matrix<double,-1,-1>::Constant(K, K,
                                std::numeric_limits<double>::quiet_NaN());
      new (&prior_Beta_loc)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(prior_Beta_loc_data__.data(),
        K, K);
      {
        std::vector<local_scalar_t__> prior_Beta_loc_flat__;
        current_statement__ = 38;
        prior_Beta_loc_flat__ = context__.vals_r("prior_Beta_loc");
        current_statement__ = 38;
        pos__ = 1;
        current_statement__ = 38;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 38;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 38;
            stan::model::assign(prior_Beta_loc, prior_Beta_loc_flat__[(pos__
              - 1)], "assigning variable prior_Beta_loc",
              stan::model::index_uni(sym2__), stan::model::index_uni(sym1__));
            current_statement__ = 38;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 39;
      stan::math::validate_non_negative_index("prior_Beta_scale", "K", K);
      current_statement__ = 40;
      stan::math::validate_non_negative_index("prior_Beta_scale", "K", K);
      current_statement__ = 41;
      context__.validate_dims("data initialization", "prior_Beta_scale",
        "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(K)});
      prior_Beta_scale_data__ = Eigen::Matrix<double,-1,-1>::Constant(K, K,
                                  std::numeric_limits<double>::quiet_NaN());
      new (&prior_Beta_scale)
        Eigen::Map<Eigen::Matrix<double,-1,-1>>(prior_Beta_scale_data__.data(),
        K, K);
      {
        std::vector<local_scalar_t__> prior_Beta_scale_flat__;
        current_statement__ = 41;
        prior_Beta_scale_flat__ = context__.vals_r("prior_Beta_scale");
        current_statement__ = 41;
        pos__ = 1;
        current_statement__ = 41;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 41;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 41;
            stan::model::assign(prior_Beta_scale,
              prior_Beta_scale_flat__[(pos__ - 1)],
              "assigning variable prior_Beta_scale",
              stan::model::index_uni(sym2__), stan::model::index_uni(sym1__));
            current_statement__ = 41;
            pos__ = (pos__ + 1);
          }
        }
      }
      current_statement__ = 42;
      context__.validate_dims("data initialization", "prior_Rho_marginal",
        "int", std::vector<size_t>{});
      prior_Rho_marginal = std::numeric_limits<int>::min();
      current_statement__ = 42;
      prior_Rho_marginal = context__.vals_i("prior_Rho_marginal")[(1 - 1)];
      current_statement__ = 42;
      stan::math::check_greater_or_equal(function__, "prior_Rho_marginal",
        prior_Rho_marginal, 1);
      current_statement__ = 43;
      stan::math::validate_non_negative_index("I", "K", K);
      current_statement__ = 44;
      stan::math::validate_non_negative_index("I", "K", K);
      current_statement__ = 45;
      I_data__ = Eigen::Matrix<double,-1,-1>::Constant(K, K,
                   std::numeric_limits<double>::quiet_NaN());
      new (&I) Eigen::Map<Eigen::Matrix<double,-1,-1>>(I_data__.data(), K, K);
      current_statement__ = 45;
      stan::model::assign(I,
        stan::math::diag_matrix(stan::math::rep_vector(1, K)),
        "assigning variable I");
      current_statement__ = 46;
      stan::math::validate_non_negative_index("Beta_raw", "K", K);
      current_statement__ = 47;
      stan::math::validate_non_negative_index("Beta_raw", "K", K);
      current_statement__ = 48;
      stan::math::validate_non_negative_index("Theta", "K", K);
      current_statement__ = 48;
      stan::math::validate_non_negative_index("Theta", "K", K);
      current_statement__ = 49;
      stan::math::validate_non_negative_index("Beta", "K", K);
      current_statement__ = 50;
      stan::math::validate_non_negative_index("Beta", "K", K);
      current_statement__ = 51;
      stan::math::validate_non_negative_index("Sigma", "K", K);
      current_statement__ = 52;
      stan::math::validate_non_negative_index("Sigma", "K", K);
      current_statement__ = 53;
      stan::math::validate_non_negative_index("Rho", "K", K);
      current_statement__ = 54;
      stan::math::validate_non_negative_index("Rho", "K", K);
      current_statement__ = 55;
      log_lik_1dim__ = std::numeric_limits<int>::min();
      current_statement__ = 55;
      log_lik_1dim__ = (T - 1);
      current_statement__ = 55;
      stan::math::validate_non_negative_index("log_lik", "T - 1",
        log_lik_1dim__);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    num_params_r__ = (K * K) + (K + ((K * (K - 1)) / 2));
  }
  inline std::string model_name() const final {
    return "model_VAR_wishart";
  }
  inline std::vector<std::string> model_compile_info() const noexcept {
    return std::vector<std::string>{"stanc_version = stanc3 v2.32.2",
             "stancflags = --allow-undefined"};
  }
  template <bool propto__, bool jacobian__, typename VecR, typename VecI,
            stan::require_vector_like_t<VecR>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline stan::scalar_type_t<VecR>
  log_prob_impl(VecR& params_r__, VecI& params_i__, std::ostream*
                pstream__ = nullptr) const {
    using T__ = stan::scalar_type_t<VecR>;
    using local_scalar_t__ = T__;
    T__ lp__(0.0);
    stan::math::accumulator<T__> lp_accum__;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    static constexpr const char* function__ =
      "model_VAR_wishart_namespace::log_prob";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<local_scalar_t__,-1,-1> Beta_raw =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      current_statement__ = 1;
      Beta_raw = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K,
                   K);
      Eigen::Matrix<local_scalar_t__,-1,-1> Theta =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      current_statement__ = 2;
      Theta = in__.template read_constrain_cov_matrix<
                Eigen::Matrix<local_scalar_t__,-1,-1>, jacobian__>(lp__, K);
      Eigen::Matrix<local_scalar_t__,-1,-1> Beta =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      current_statement__ = 3;
      stan::model::assign(Beta,
        stan::math::add(stan::math::elt_multiply(Beta_raw, prior_Beta_scale),
          prior_Beta_loc), "assigning variable Beta");
      Eigen::Matrix<local_scalar_t__,-1,-1> Sigma =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      current_statement__ = 4;
      stan::model::assign(Sigma, stan::math::inverse_spd(Theta),
        "assigning variable Sigma");
      Eigen::Matrix<local_scalar_t__,-1,-1> Rho =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      {
        current_statement__ = 15;
        for (int i = 1; i <= K; ++i) {
          current_statement__ = 13;
          for (int j = 1; j <= K; ++j) {
            current_statement__ = 11;
            if (stan::math::logical_neq(i, j)) {
              current_statement__ = 9;
              stan::model::assign(Rho,
                (-stan::model::rvalue(Theta, "Theta",
                    stan::model::index_uni(i), stan::model::index_uni(j)) /
                stan::math::sqrt(
                  (stan::model::rvalue(Theta, "Theta",
                     stan::model::index_uni(i), stan::model::index_uni(i)) *
                  stan::model::rvalue(Theta, "Theta",
                    stan::model::index_uni(j), stan::model::index_uni(j))))),
                "assigning variable Rho", stan::model::index_uni(i),
                stan::model::index_uni(j));
            } else {
              current_statement__ = 7;
              stan::model::assign(Rho, 0, "assigning variable Rho",
                stan::model::index_uni(i), stan::model::index_uni(j));
            }
          }
        }
      }
      {
        current_statement__ = 23;
        lp_accum__.add(stan::math::std_normal_lpdf<false>(
                         stan::math::to_vector(Beta_raw)));
        current_statement__ = 24;
        lp_accum__.add(stan::math::inv_wishart_lpdf<false>(Theta, ((((1 /
                         prior_Rho_marginal) - 1) + K) - 1), I));
        {
          current_statement__ = 29;
          for (int t = 2; t <= T; ++t) {
            current_statement__ = 25;
            stan::math::validate_non_negative_index("mu", "K", K);
            Eigen::Matrix<local_scalar_t__,-1,1> mu =
              Eigen::Matrix<local_scalar_t__,-1,1>::Constant(K, DUMMY_VAR__);
            current_statement__ = 26;
            stan::model::assign(mu,
              stan::math::multiply(Beta,
                stan::model::rvalue(Y, "Y", stan::model::index_uni((t - 1)),
                  stan::model::index_omni())), "assigning variable mu");
            current_statement__ = 27;
            lp_accum__.add(stan::math::multi_normal_lpdf<false>(
                             stan::model::rvalue(Y, "Y",
                               stan::model::index_uni(t),
                               stan::model::index_omni()), mu, Sigma));
          }
        }
      }
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
    lp_accum__.add(lp__);
    return lp_accum__.sum();
  }
  template <typename RNG, typename VecR, typename VecI, typename VecVar,
            stan::require_vector_like_vt<std::is_floating_point,
            VecR>* = nullptr, stan::require_vector_like_vt<std::is_integral,
            VecI>* = nullptr, stan::require_vector_vt<std::is_floating_point,
            VecVar>* = nullptr>
  inline void
  write_array_impl(RNG& base_rng__, VecR& params_r__, VecI& params_i__,
                   VecVar& vars__, const bool
                   emit_transformed_parameters__ = true, const bool
                   emit_generated_quantities__ = true, std::ostream*
                   pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    static constexpr bool propto__ = true;
    // suppress unused var warning
    (void) propto__;
    double lp__ = 0.0;
    // suppress unused var warning
    (void) lp__;
    int current_statement__ = 0;
    stan::math::accumulator<double> lp_accum__;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    constexpr bool jacobian__ = false;
    static constexpr const char* function__ =
      "model_VAR_wishart_namespace::write_array";
    // suppress unused var warning
    (void) function__;
    try {
      Eigen::Matrix<double,-1,-1> Beta_raw =
        Eigen::Matrix<double,-1,-1>::Constant(K, K,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 1;
      Beta_raw = in__.template read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K,
                   K);
      Eigen::Matrix<double,-1,-1> Theta =
        Eigen::Matrix<double,-1,-1>::Constant(K, K,
          std::numeric_limits<double>::quiet_NaN());
      current_statement__ = 2;
      Theta = in__.template read_constrain_cov_matrix<
                Eigen::Matrix<local_scalar_t__,-1,-1>, jacobian__>(lp__, K);
      Eigen::Matrix<double,-1,-1> Beta =
        Eigen::Matrix<double,-1,-1>::Constant(K, K,
          std::numeric_limits<double>::quiet_NaN());
      Eigen::Matrix<double,-1,-1> Sigma =
        Eigen::Matrix<double,-1,-1>::Constant(K, K,
          std::numeric_limits<double>::quiet_NaN());
      Eigen::Matrix<double,-1,-1> Rho =
        Eigen::Matrix<double,-1,-1>::Constant(K, K,
          std::numeric_limits<double>::quiet_NaN());
      out__.write(Beta_raw);
      out__.write(Theta);
      if (stan::math::logical_negation(
            (stan::math::primitive_value(emit_transformed_parameters__) ||
            stan::math::primitive_value(emit_generated_quantities__)))) {
        return ;
      }
      current_statement__ = 3;
      stan::model::assign(Beta,
        stan::math::add(stan::math::elt_multiply(Beta_raw, prior_Beta_scale),
          prior_Beta_loc), "assigning variable Beta");
      current_statement__ = 4;
      stan::model::assign(Sigma, stan::math::inverse_spd(Theta),
        "assigning variable Sigma");
      {
        current_statement__ = 15;
        for (int i = 1; i <= K; ++i) {
          current_statement__ = 13;
          for (int j = 1; j <= K; ++j) {
            current_statement__ = 11;
            if (stan::math::logical_neq(i, j)) {
              current_statement__ = 9;
              stan::model::assign(Rho,
                (-stan::model::rvalue(Theta, "Theta",
                    stan::model::index_uni(i), stan::model::index_uni(j)) /
                stan::math::sqrt(
                  (stan::model::rvalue(Theta, "Theta",
                     stan::model::index_uni(i), stan::model::index_uni(i)) *
                  stan::model::rvalue(Theta, "Theta",
                    stan::model::index_uni(j), stan::model::index_uni(j))))),
                "assigning variable Rho", stan::model::index_uni(i),
                stan::model::index_uni(j));
            } else {
              current_statement__ = 7;
              stan::model::assign(Rho, 0, "assigning variable Rho",
                stan::model::index_uni(i), stan::model::index_uni(j));
            }
          }
        }
      }
      if (emit_transformed_parameters__) {
        out__.write(Beta);
        out__.write(Sigma);
        out__.write(Rho);
      }
      if (stan::math::logical_negation(emit_generated_quantities__)) {
        return ;
      }
      Eigen::Matrix<double,-1,1> log_lik =
        Eigen::Matrix<double,-1,1>::Constant(log_lik_1dim__,
          std::numeric_limits<double>::quiet_NaN());
      {
        current_statement__ = 21;
        for (int t = 2; t <= T; ++t) {
          current_statement__ = 17;
          stan::math::validate_non_negative_index("mu", "K", K);
          Eigen::Matrix<double,-1,1> mu =
            Eigen::Matrix<double,-1,1>::Constant(K,
              std::numeric_limits<double>::quiet_NaN());
          current_statement__ = 18;
          stan::model::assign(mu,
            stan::math::multiply(Beta,
              stan::model::rvalue(Y, "Y", stan::model::index_uni((t - 1)),
                stan::model::index_omni())), "assigning variable mu");
          current_statement__ = 19;
          stan::model::assign(log_lik,
            stan::math::multi_normal_lpdf<false>(
              stan::model::rvalue(Y, "Y", stan::model::index_uni(t),
                stan::model::index_omni()), mu, Sigma),
            "assigning variable log_lik", stan::model::index_uni((t - 1)));
        }
      }
      out__.write(log_lik);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, typename VecI,
            stan::require_vector_t<VecVar>* = nullptr,
            stan::require_vector_like_vt<std::is_integral, VecI>* = nullptr>
  inline void
  unconstrain_array_impl(const VecVar& params_r__, const VecI& params_i__,
                         VecVar& vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::deserializer<local_scalar_t__> in__(params_r__, params_i__);
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,-1> Beta_raw =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      current_statement__ = 1;
      stan::model::assign(Beta_raw,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K, K),
        "assigning variable Beta_raw");
      out__.write(Beta_raw);
      Eigen::Matrix<local_scalar_t__,-1,-1> Theta =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      current_statement__ = 2;
      stan::model::assign(Theta,
        in__.read<Eigen::Matrix<local_scalar_t__,-1,-1>>(K, K),
        "assigning variable Theta");
      out__.write_free_cov_matrix(Theta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  template <typename VecVar, stan::require_vector_t<VecVar>* = nullptr>
  inline void
  transform_inits_impl(const stan::io::var_context& context__, VecVar&
                       vars__, std::ostream* pstream__ = nullptr) const {
    using local_scalar_t__ = double;
    stan::io::serializer<local_scalar_t__> out__(vars__);
    int current_statement__ = 0;
    local_scalar_t__ DUMMY_VAR__(std::numeric_limits<double>::quiet_NaN());
    // suppress unused var warning
    (void) DUMMY_VAR__;
    try {
      current_statement__ = 1;
      context__.validate_dims("parameter initialization", "Beta_raw",
        "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(K)});
      current_statement__ = 2;
      context__.validate_dims("parameter initialization", "Theta", "double",
        std::vector<size_t>{static_cast<size_t>(K), static_cast<size_t>(K)});
      int pos__ = std::numeric_limits<int>::min();
      pos__ = 1;
      Eigen::Matrix<local_scalar_t__,-1,-1> Beta_raw =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> Beta_raw_flat__;
        current_statement__ = 1;
        Beta_raw_flat__ = context__.vals_r("Beta_raw");
        current_statement__ = 1;
        pos__ = 1;
        current_statement__ = 1;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 1;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 1;
            stan::model::assign(Beta_raw, Beta_raw_flat__[(pos__ - 1)],
              "assigning variable Beta_raw", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 1;
            pos__ = (pos__ + 1);
          }
        }
      }
      out__.write(Beta_raw);
      Eigen::Matrix<local_scalar_t__,-1,-1> Theta =
        Eigen::Matrix<local_scalar_t__,-1,-1>::Constant(K, K, DUMMY_VAR__);
      {
        std::vector<local_scalar_t__> Theta_flat__;
        current_statement__ = 2;
        Theta_flat__ = context__.vals_r("Theta");
        current_statement__ = 2;
        pos__ = 1;
        current_statement__ = 2;
        for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
          current_statement__ = 2;
          for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
            current_statement__ = 2;
            stan::model::assign(Theta, Theta_flat__[(pos__ - 1)],
              "assigning variable Theta", stan::model::index_uni(sym2__),
              stan::model::index_uni(sym1__));
            current_statement__ = 2;
            pos__ = (pos__ + 1);
          }
        }
      }
      out__.write_free_cov_matrix(Theta);
    } catch (const std::exception& e) {
      stan::lang::rethrow_located(e, locations_array__[current_statement__]);
    }
  }
  inline void
  get_param_names(std::vector<std::string>& names__, const bool
                  emit_transformed_parameters__ = true, const bool
                  emit_generated_quantities__ = true) const {
    names__ = std::vector<std::string>{"Beta_raw", "Theta"};
    if (emit_transformed_parameters__) {
      std::vector<std::string> temp{"Beta", "Sigma", "Rho"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::string> temp{"log_lik"};
      names__.reserve(names__.size() + temp.size());
      names__.insert(names__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  get_dims(std::vector<std::vector<size_t>>& dimss__, const bool
           emit_transformed_parameters__ = true, const bool
           emit_generated_quantities__ = true) const {
    dimss__ = std::vector<std::vector<size_t>>{std::vector<size_t>{static_cast<
                                                                    size_t>(K),
                                                 static_cast<size_t>(K)},
                std::vector<size_t>{static_cast<size_t>(K),
                  static_cast<size_t>(K)}};
    if (emit_transformed_parameters__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(K),
               static_cast<size_t>(K)},
             std::vector<size_t>{static_cast<size_t>(K),
               static_cast<size_t>(K)},
             std::vector<size_t>{static_cast<size_t>(K),
               static_cast<size_t>(K)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
    if (emit_generated_quantities__) {
      std::vector<std::vector<size_t>>
        temp{std::vector<size_t>{static_cast<size_t>(log_lik_1dim__)}};
      dimss__.reserve(dimss__.size() + temp.size());
      dimss__.insert(dimss__.end(), temp.begin(), temp.end());
    }
  }
  inline void
  constrained_param_names(std::vector<std::string>& param_names__, bool
                          emit_transformed_parameters__ = true, bool
                          emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "Beta_raw" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "Theta" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "Beta" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "Sigma" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "Rho" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
    }
  }
  inline void
  unconstrained_param_names(std::vector<std::string>& param_names__, bool
                            emit_transformed_parameters__ = true, bool
                            emit_generated_quantities__ = true) const final {
    for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
      for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
        param_names__.emplace_back(std::string() + "Beta_raw" + '.' +
          std::to_string(sym2__) + '.' + std::to_string(sym1__));
      }
    }
    for (int sym1__ = 1; sym1__ <= (K + ((K * (K - 1)) / 2)); ++sym1__) {
      param_names__.emplace_back(std::string() + "Theta" + '.' +
        std::to_string(sym1__));
    }
    if (emit_transformed_parameters__) {
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "Beta" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "Sigma" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
      for (int sym1__ = 1; sym1__ <= K; ++sym1__) {
        for (int sym2__ = 1; sym2__ <= K; ++sym2__) {
          param_names__.emplace_back(std::string() + "Rho" + '.' +
            std::to_string(sym2__) + '.' + std::to_string(sym1__));
        }
      }
    }
    if (emit_generated_quantities__) {
      for (int sym1__ = 1; sym1__ <= log_lik_1dim__; ++sym1__) {
        param_names__.emplace_back(std::string() + "log_lik" + '.' +
          std::to_string(sym1__));
      }
    }
  }
  inline std::string get_constrained_sizedtypes() const {
    return std::string("[{\"name\":\"Beta_raw\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"parameters\"},{\"name\":\"Theta\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"parameters\"},{\"name\":\"Beta\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"transformed_parameters\"},{\"name\":\"Sigma\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"transformed_parameters\"},{\"name\":\"Rho\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(log_lik_1dim__) + "},\"block\":\"generated_quantities\"}]");
  }
  inline std::string get_unconstrained_sizedtypes() const {
    return std::string("[{\"name\":\"Beta_raw\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"parameters\"},{\"name\":\"Theta\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string((K + ((K * (K - 1)) /2))) + "},\"block\":\"parameters\"},{\"name\":\"Beta\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"transformed_parameters\"},{\"name\":\"Sigma\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"transformed_parameters\"},{\"name\":\"Rho\",\"type\":{\"name\":\"matrix\",\"rows\":" + std::to_string(K) + ",\"cols\":" + std::to_string(K) + "},\"block\":\"transformed_parameters\"},{\"name\":\"log_lik\",\"type\":{\"name\":\"vector\",\"length\":" + std::to_string(log_lik_1dim__) + "},\"block\":\"generated_quantities\"}]");
  }
  // Begin method overload boilerplate
  template <typename RNG> inline void
  write_array(RNG& base_rng, Eigen::Matrix<double,-1,1>& params_r,
              Eigen::Matrix<double,-1,1>& vars, const bool
              emit_transformed_parameters = true, const bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((K * K) + (K * K));
    const size_t num_transformed = emit_transformed_parameters * ((((K * K) +
      (K * K)) + (K * K)));
    const size_t num_gen_quantities = emit_generated_quantities *
      (log_lik_1dim__);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    std::vector<int> params_i;
    vars = Eigen::Matrix<double,-1,1>::Constant(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <typename RNG> inline void
  write_array(RNG& base_rng, std::vector<double>& params_r, std::vector<int>&
              params_i, std::vector<double>& vars, bool
              emit_transformed_parameters = true, bool
              emit_generated_quantities = true, std::ostream*
              pstream = nullptr) const {
    const size_t num_params__ = ((K * K) + (K * K));
    const size_t num_transformed = emit_transformed_parameters * ((((K * K) +
      (K * K)) + (K * K)));
    const size_t num_gen_quantities = emit_generated_quantities *
      (log_lik_1dim__);
    const size_t num_to_write = num_params__ + num_transformed +
      num_gen_quantities;
    vars = std::vector<double>(num_to_write,
             std::numeric_limits<double>::quiet_NaN());
    write_array_impl(base_rng, params_r, params_i, vars,
      emit_transformed_parameters, emit_generated_quantities, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(Eigen::Matrix<T_,-1,1>& params_r, std::ostream* pstream = nullptr) const {
    Eigen::Matrix<int,-1,1> params_i;
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  template <bool propto__, bool jacobian__, typename T_> inline T_
  log_prob(std::vector<T_>& params_r, std::vector<int>& params_i,
           std::ostream* pstream = nullptr) const {
    return log_prob_impl<propto__, jacobian__>(params_r, params_i, pstream);
  }
  inline void
  transform_inits(const stan::io::var_context& context,
                  Eigen::Matrix<double,-1,1>& params_r, std::ostream*
                  pstream = nullptr) const final {
    std::vector<double> params_r_vec(params_r.size());
    std::vector<int> params_i;
    transform_inits(context, params_i, params_r_vec, pstream);
    params_r = Eigen::Map<Eigen::Matrix<double,-1,1>>(params_r_vec.data(),
                 params_r_vec.size());
  }
  inline void
  transform_inits(const stan::io::var_context& context, std::vector<int>&
                  params_i, std::vector<double>& vars, std::ostream*
                  pstream__ = nullptr) const {
    vars.resize(num_params_r__);
    transform_inits_impl(context, vars, pstream__);
  }
  inline void
  unconstrain_array(const std::vector<double>& params_constrained,
                    std::vector<double>& params_unconstrained, std::ostream*
                    pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = std::vector<double>(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
  inline void
  unconstrain_array(const Eigen::Matrix<double,-1,1>& params_constrained,
                    Eigen::Matrix<double,-1,1>& params_unconstrained,
                    std::ostream* pstream = nullptr) const {
    const std::vector<int> params_i;
    params_unconstrained = Eigen::Matrix<double,-1,1>::Constant(num_params_r__,
                             std::numeric_limits<double>::quiet_NaN());
    unconstrain_array_impl(params_constrained, params_i,
      params_unconstrained, pstream);
  }
};
}
using stan_model = model_VAR_wishart_namespace::model_VAR_wishart;
#ifndef USING_R
// Boilerplate
stan::model::model_base&
new_model(stan::io::var_context& data_context, unsigned int seed,
          std::ostream* msg_stream) {
  stan_model* m = new stan_model(data_context, seed, msg_stream);
  return *m;
}
stan::math::profile_map& get_stan_profile_data() {
  return model_VAR_wishart_namespace::profiles__;
}
#endif
#endif
