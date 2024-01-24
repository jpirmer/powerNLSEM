### powerNLSEM 0.1.0

* add Wald-test/two sided z-test to significance decisions which can be modeled (additionally to one sided z-tests)
* add factor score regression to model estimation techniques
* add unconstrained product indicator approach to model estimation techniques
* add estimates and ses in saved objects which are more easily reanalyzed
* add possibility to change the relation to `n` (default is `sqrt(n)`, which is the theoretical justified and also performed best in simulation studies)
* Non-convergences and possible faulty model fits are now marked and taken out of the estimation of power

<!--- check that these are actually implemented!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! -->

### powerNLSEM 0.0.2

* add plot and reanalyzation to `powerNLSEM` objects
* add both probit and logit regression to be used in the estimation of power


### powerNLSEM 0.0.1

* Initial commit of the `powerNLSEM` package with first working functions

