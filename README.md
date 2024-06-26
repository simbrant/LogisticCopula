 An implementation of a method of extending a logistic regression
    model beyond linear effects of the covariates. The extension in is
    constructed by first equating the logistic regression model to a naive Bayes
    model where all the margins are specfied to follow natural exponential 
    distributions conditional on Y, that is, a model for Y given X that is
    specified through the distribution of X given Y, where the columns of X are
    assumed to be mutually independent conditional on Y. Subsequently, the
    model is expanded by adding vine - copulas to relax the assumption of
    mutual independence, where pair-copulas are added in a stage-wise, forward
    selection manner. Some heuristics are employed during the process of
    selecting edges, as well as the families of pair-copula models. After each
    component is added, the parameters are updated by a (smaller) number of
    gradient steps to maximise the likelihood. When the algorithm has stopped
    adding edges, based the criterion that a new edge should improve the
    likelihood more than k times the number new parameters, the parameters are
    updated with a larger number of gradient steps, or until convergence. 
