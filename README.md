# FREEtree

FREEtree, a tree-based method for high dimensional longitudinal data with correlated features. FREEtree deals with longitudinal data by using a  piecewise random effect model. It also exploits the network structure of the features, by first clustering them using Weighted Gene Coexpression Network Analysis (WGCNA). Then it conducts a screening step within each cluster of features and a selecting step among the surviving features, which provides a relatively unbiased way to select features. By using dominant principle components as regression variables at each leaf and the original features as splitting variables at splitting nodes, FREEtree maintains its interpretability and improves its computational efficiency. The simulation results show that FREEtree has improvements over other tree-based methods in terms of prediction accuracy, feature selection accuracy as well as the ability to recover the underlying correct structure. 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

### Prerequisites

What things you need to install the software and how to install them

```
Give examples
```

### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc

Source of template for this readme file: https://gist.github.com/PurpleBooth/109311bb0361f32d87a2
