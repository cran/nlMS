\name{model_selection}
\alias{model_selection}
\title{Model Selection (MS)}
\usage{
model_selection(data=NULL, name_dir=NULL,
                name_results=NULL)
}
\description{
This function fits the input data to nine regression models. It gives AIC (Akaike Information Index) and RSE (Residual Standard Error) as goodness of fit measures. The models implemented are: "linear", "quadratic", "cubic", "logistic", "exponential", "power", "monod", "haldane", "logit".
For the function gnls() from the 'nlme' package is taken for all the models.
Furthermore, the shapiro.test() function from the 'stats' package is used to test residuals normality.
The function has been though to run over several values. This is, when 'data' is composed by several dependent variables Y.
}
\arguments{
\item{data}{Matrix with at least two columns. Fisrt column is taken as the X, ie the independent variable. Next columns are recursevelly taken as the Y variable (dependent). Note: in the current version X must be positive.}
\item{name_dir}{Directory to save function results. Default is a temporary directory set by tempdir().}
\item{name_results}{Name of the file for the global results. Default is "nlMS_results". Note: by default, the extension of the output file is .csv. You do not need to specify it.}
}
\value{There are two kinds of outputs. First one is a .csv file, that takes its name from the argument name_results, it consists of a table with the value of AIC, RSE and R.squared (only for linear models) for each tested model. As well a fourth row with the Saphiro.Wilk normality test result applied to the residuals. If you want to use AIC, then the resdiuals must follow a normal distribution.
In case input Y data consists of more than one column, the model performs several runs always taking the same variable X as independent values.}
\author{Carme Font <carme.font.moragon@gmail.com>}
\references{
The package has been tested in

Romero, F., Sabater, S., Font, C., Balcazar, J. L., Acuna, V. (2019). Desiccation events change the microbial response to gradients of wastewater effluent pollution. Water research, 151, 371-380. doi: 10.1016/j.watres.2018.12.028

Acuna, V., Casellas, M., Font, C., Romero, F., Sabater, S. (2019). Nutrient attenuation dynamics in effluent dominated watercourses. Water Research. doi: 10.1016/j.watres.2019.05.093

Sabater-Liesa, L., Montemurro, N., Font, C., Ginebreda, A., Gonzalez-Trujillo, J. D., Mingorance, N., ... Barcelo, D. (2019). The response patterns of stream biofilms to urban sewage change with exposure time and dilution. Science of The Total Environment, 674, 401-411. doi: 10.1016/j.scitotenv.2019.04.178
}
\examples{
X<-rnorm(100)
X<-X-min(X)+0.001
Y<-rnorm(100)
model_selection(data=cbind(X, Y))
}
