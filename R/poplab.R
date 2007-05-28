"age.profile" <-
function (ped, sex, year) 
{
    if (!(sex %in% c(1, 2))) 
        stop("Gender not recognized")
    ped <- ped[(ped[, "yod"] == 0 | ped[, "yod"] > year) & ped[, 
        "yob"] <= year & ped[, "sex"] == sex, ]
    ped <- cbind(ped, age = year - ped[, "yob"])
    runinprof <- cbind(as.integer(rownames(as.matrix(table(ped[, 
        "age"])))), as.matrix(table(ped[, "age"])))
    colnames(runinprof) <- c("age", "counts")
    return(runinprof)
}
"assigncancer" <-
function (ped, inc, risk, yr, sex, type, fam.rel) 
{
    if (!(sex %in% c(1, 2))) 
        stop("Gender not recognized")
    if (!(fam.rel %in% c("s", "p"))) 
        stop("Familial relationship for disease aggregation not recognized")
    condped <- ped[ped[, "yod"] == 0 & ped[, "yoi"] == 0 & ped[, 
        "sex"] == sex, ]
    condped <- cbind(condped, age = (yr - condped[, "yob"]))
    condped <- cbind(condped, rate = inc[, "rate"][match(condped[, 
        "age"], inc[, "age"])])
    condped <- condped[!is.na(condped[, "rate"]), ]
    if (fam.rel == "s") {
        condped <- condped[condped[, "m"] != 0, ]
        ped <- ped[order(ped[, "yoi"]), ]
        condped <- cbind(condped, yoi.rel = ped[ped[, "yoi"] != 
            0, "yoi"][match(condped[, "m"], ped[ped[, "yoi"] != 
            0, "m"])], age.i.rel = ped[ped[, "yoi"] != 0, "yob"][match(condped[, 
            "m"], ped[ped[, "yoi"] != 0, "m"])])
    }
    else {
        if (sex == 1) {
            condped <- cbind(condped, yoi.rel = ped[, "yoi"][match(condped[, 
                "f"], ped[, "ID"])], age.i.rel = ped[, "yob"][match(condped[, 
                "f"], ped[, "ID"])])
        }
        else {
            condped <- cbind(condped, yoi.rel = ped[, "yoi"][match(condped[, 
                "m"], ped[, "ID"])], age.i.rel = ped[, "yob"][match(condped[, 
                "m"], ped[, "ID"])])
        }
    }
    condped[is.na(condped[, "yoi.rel"]), "yoi.rel"] <- 0
    condped[, "age.i.rel"] <- ifelse(condped[, "yoi.rel"] == 
        0, 0, condped[, "yoi.rel"] - condped[, "age.i.rel"])
    condped[(condped[, "yoi.rel"] != 0), "rate"] <- switch(type, 
        rr = condped[(condped[, "yoi.rel"] != 0), "rate"] * risk, 
        or = risk * condped[(condped[, "yoi.rel"] != 0), "rate"]/(condped[(condped[, 
            "yoi.rel"] != 0), "rate"] * (risk - 1) + 1), agesprr = condped[condped[, 
            "yoi.rel"] != 0, "rate"] * risk[cut(condped[condped[, 
            "yoi.rel"] != 0, "age.i.rel"], c(0, risk[, "age"]), 
            right = FALSE, include.lowest = TRUE), "val"], agespor = risk[cut(condped[condped[, 
            "yoi.p"] != 0, "age.i.rel"], c(0, risk[, "age"]), 
            right = FALSE, include.lowest = TRUE), "val"] * condped[(condped[, 
            "yoi.rel"] != 0), "rate"]/(condped[(condped[, "yoi.rel"] != 
            0), "rate"] * (risk[cut(condped[condped[, "yoi.rel"] != 
            0, "age.i.rel"], c(0, risk[, "age"]), right = FALSE, 
            include.lowest = TRUE), "val"] - 1) + 1), )
    if (nrow(condped[condped[, "rate"] > 1, , drop = FALSE]) > 
        0) {
        warning("One (or more) probability of disease incidence greater than 1 \n")
        condped[condped[, "rate"] > 1, "rate"] <- 1
    }
    condped <- cbind(condped, cancer = rbinom(nrow(condped), 
        1, condped[, "rate"]))
    ped[ped[, "ID"] %in% condped[condped[, "cancer"] == 1, "ID"], 
        "yoi"] <- yr
    return(ped)
}
"assigndeath" <-
function (ped, tempmort, risk, yr) 
{
    condped <- cbind(ped[ped[, "yod"] == 0, , drop = FALSE], 
        age = (yr - ped[ped[, "yod"] == 0, "yob"]))
    condped <- rbind(cbind(condped[condped[, "sex"] == 2, ], 
        rate = tempmort[, "femrate"][match(condped[condped[, 
            "sex"] == 2, "age"], tempmort[, "age"])]), cbind(condped[condped[, 
        "sex"] == 1, ], rate = tempmort[, "malerate"][match(condped[condped[, 
        "sex"] == 1, "age"], tempmort[, "age"])]))
    condped <- condped[!is.na(condped[, "rate"]) & condped[, 
        "rate"] != 0, ]
    condped[condped[, "yoi"] != 0, "rate"] <- condped[condped[, 
        "yoi"] != 0, "rate"] * risk
    if (nrow(condped[condped[, "rate"] > 1, , drop = FALSE]) > 
        0) {
        warning("One (or more) probability of death greater than 1 \n")
        condped[condped[, "rate"] > 1, "rate"] <- 1
    }
    condped <- cbind(condped, mortch = rbinom(nrow(condped), 
        1, condped[, "rate"]))
    ped[(ped[, "ID"] %in% condped[condped[, "mortch"] == 1, "ID"]) | 
        (ped[, "ID"] %in% condped[condped[, "age"] >= 100, "ID"]), 
        "yod"] <- yr
    return(ped)
}
"create.baseline.complete" <-
function (baseyear, healthy, risk, sex.a, mortratio, base.scale, 
    runintime, d.mod, fam.rel, print.option, population.fem, 
    population.male, mortality.fem, mortality.male, fertility, 
    incidence, seed, folder) 
{
    if (missing(healthy)) {
        healthy <- 1
    }
    if (missing(risk)) {
        risk <- 1
    }
    if (missing(sex.a)) {
        sex.a <- 2
    }
    if (missing(mortratio)) {
        mortratio <- 1
    }
    if (missing(base.scale)) {
        base.scale <- 500
    }
    if (missing(runintime)) {
        runintime <- 100
    }
    if (missing(print.option)) {
        print.option <- FALSE
    }
    if (missing(seed)) {
        seed <- NULL
    }
    if ((!is.null(seed)) & (!is.character(seed))) 
        set.seed(seed)
    if (base.scale <= 0) 
        stop("Something wrong with parameter base.scale. See help(create.baseline.complete)")
    if (runintime <= 0) 
        stop("Something wrong with parameter runintime. See help(create.baseline.complete)")
    baseyr <- baseyear - runintime
    fert <- readin(f.name = fertility, folder)
    if (!(baseyear %in% as.integer(colnames(fert)[-1]))) 
        stop(paste("No fertility rates for the year", baseyear, 
            "in data files"))
    childbc <- fert[, c(1, which(colnames(fert) == baseyear, 
        arr.ind = TRUE))]
    childbc <- childbc[childbc[, 2] != 0, ]
    colnames(childbc) <- c("age", "rate")
    femmort <- readin(f.name = mortality.fem, folder)
    malemort <- readin(f.name = mortality.male, folder)
    if ((!(baseyear %in% as.integer(colnames(femmort)[-1]))) | 
        (!(baseyear %in% as.integer(colnames(malemort)[-1])))) 
        stop(paste("No mortality rates either for males or females, for the year", 
            baseyear, "in data files"))
    if (nrow(femmort) != nrow(malemort)) 
        stop("Different age structure in the male and female mortality files. Organize files and try again.")
    tempmort <- cbind(femmort[, c(1, which(colnames(femmort) == 
        baseyear, arr.ind = TRUE))], malemort[, which(colnames(malemort) == 
        baseyear, arr.ind = TRUE)])
    colnames(tempmort) <- c("age", "femrate", "malerate")
    fempop <- readin(f.name = population.fem, folder)
    if (!(baseyear %in% as.integer(colnames(fempop)[-1]))) 
        stop(paste("No female population counts for the year", 
            baseyear, "in data files"))
    fempop <- fempop[, c(1, which(colnames(fempop) == baseyear, 
        arr.ind = TRUE))]
    colnames(fempop) <- c("age", "counts")
    malepop <- readin(f.name = population.male, folder)
    if (!(baseyear %in% as.integer(colnames(malepop)[-1]))) 
        stop(paste("No male population counts for the year", 
            baseyear, "in data files"))
    malepop <- malepop[, c(1, which(colnames(malepop) == baseyear, 
        arr.ind = TRUE))]
    colnames(malepop) <- c("age", "counts")
    if (healthy == 0) {
        inc <- readin(f.name = incidence, folder)
        if (!(baseyear %in% as.integer(colnames(inc)[-1]))) 
            stop(paste("No incidence rates for the year", baseyear, 
                "in data files"))
        inc <- inc[, c(1, which(colnames(inc) == baseyear, arr.ind = TRUE))]
        inc <- inc[inc[, 2] != 0, , drop = FALSE]
        colnames(inc) <- c("age", "rate")
        if (length(risk) > 1) {
            risk <- matrix(data = risk, byrow = FALSE, ncol = 2)
            colnames(risk) <- c("age", "val")
            if (!(d.mod %in% c("agesprr", "agespor"))) 
                stop(paste("Invalid risk value specification for the chosen disease model"))
        }
        else {
            if (!(d.mod %in% c("rr", "or"))) 
                stop(paste("Invalid risk value specification for the chosen disease model"))
        }
        if (!(d.mod %in% c("agespor", "agesprr", "or", "rr"))) 
            stop(paste("Invalid disease model specification"))
    }
    baseyrtot.fem <- sum(fempop[, colnames(fempop) == "counts"])
    femprof <- cbind(fempop[, "age"], fempop[, colnames(fempop) == 
        "counts"]/baseyrtot.fem * 100, (baseyr - fempop[, "age"]))
    colnames(femprof) <- c("age", "counts", "yob")
    femprof[, "counts"] <- round(femprof[, "counts"] * base.scale)
    baseyrtot.male <- sum(malepop[, colnames(malepop) == "counts"])
    maleprof <- cbind(malepop[, "age"], malepop[, colnames(malepop) == 
        "counts"]/baseyrtot.male * 100, (baseyr - malepop[, "age"]))
    colnames(maleprof) <- c("age", "counts", "yob")
    maleprof[, "counts"] <- round(maleprof[, "counts"] * base.scale * 
        (sum(malepop[, "counts"])/sum(fempop[, "counts"])))
    femprof <- femprof[(femprof[, "age"] != 0) & (femprof[, "counts"] > 
        0), ]
    maleprof <- maleprof[(maleprof[, "age"] != 0) & (maleprof[, 
        "counts"] > 0), ]
    s.pedigree <- NULL
    counter <- 0
    s.pedigree <- rbind(cbind(c(mapply(seq, cumsum(maleprof[, 
        "counts"]) - maleprof[, "counts"] + 1, cumsum(maleprof[, 
        "counts"])), recursive = TRUE), c(mapply(rep, maleprof[, 
        "yob"], maleprof[, "counts"]), recursive = TRUE), rep(1, 
        length = sum(maleprof[, "counts"])), c(mapply(rep, 0, 
        maleprof[, "counts"]), recursive = TRUE)), cbind(c(mapply(seq, 
        cumsum(femprof[, "counts"]) - femprof[, "counts"] + 1, 
        cumsum(femprof[, "counts"])), recursive = TRUE), c(mapply(rep, 
        femprof[, "yob"], femprof[, "counts"]), recursive = TRUE), 
        rep(2, length = sum(femprof[, "counts"])), c(mapply(rep, 
            0, femprof[, "counts"]), recursive = TRUE)))
    s.pedigree <- cbind(s.pedigree, rep(0, length = nrow(s.pedigree)), 
        rep(0, length = nrow(s.pedigree)), rep(0, length = nrow(s.pedigree)))
    colnames(s.pedigree) <- c("ID", "yob", "sex", "m", "f", "yod", 
        "yoi")
    correct <- max(s.pedigree[s.pedigree[, "sex"] == 1, "ID"])
    s.pedigree[s.pedigree[, "sex"] == 2, "ID"] <- s.pedigree[s.pedigree[, 
        "sex"] == 2, "ID"] + correct
    counter <- max(s.pedigree[, "ID"])
    cat("Start simulating baseline population of year", baseyear, 
        "\n")
    for (z in 1:(runintime + 1)) {
        s.pedigree <- givebirth(ped = s.pedigree, childbc = childbc, 
            nr = counter, yr = baseyr)
        counter <- max(s.pedigree[, "ID"])
        if (z <= runintime) 
            s.pedigree <- assigndeath(ped = s.pedigree, tempmort = tempmort, 
                risk = mortratio, yr = baseyr)
        if ((healthy == 0) && (z <= runintime)) 
            s.pedigree <- assigncancer(ped = s.pedigree, inc = inc, 
                risk = risk, yr = baseyr, sex = sex.a, type = d.mod, 
                fam.rel = fam.rel)
        baseyr <- baseyr + 1
    }
    s.ped <- s.pedigree[s.pedigree[, "yod"] == 0, ]
    pop <- rbind(decide(s.ped, sex = 1, malepop, baseyear), decide(s.ped, 
        sex = 2, fempop, baseyear))
    cat("Finished simulating the baseline population of year", 
        baseyear, "\n")
    notaround <- retrieve.parents(pop, s.pedigree)
    pop <- rbind(pop, notaround)
    class(pop) <- "poplab"
    if (print.option) 
        print.poplab(pop, "base", baseyear, folder)
    return(pop)
}
"decide" <-
function (ped, sex, pop, baseyr) 
{
    if (!(sex %in% c(1, 2))) 
        stop("Gender not recognized")
    ped <- ped[ped[, "sex"] == sex, ]
    pedpop <- age.profile(ped, sex, baseyr)
    pedpop <- pedpop[pedpop[, "age"] %in% pop[, "age"], ]
    pop <- pop[pop[, "age"] %in% pedpop[, "age"], ]
    this.age <- pop[pop[, "counts"] == max(pop[, "counts"]), 
        "age"]
    pop <- cbind(pop, ratio = pop[, "counts"]/pop[pop[, "counts"] == 
        max(pop[, "counts"]), "counts"])
    colnames(pop) <- c("age", "counts", "ratio")
    pedpop <- cbind(pedpop, ratio = pedpop[, "counts"]/pedpop[pedpop[, 
        "age"] == this.age, "counts"])
    colnames(pedpop) <- c("age", "counts", "ratio")
    old.value.this <- pedpop[pedpop[, "age"] == this.age, "counts"]
    compare <- cbind(pedpop[, "age"], pedpop[, "ratio"], pop[, 
        "ratio"])
    compare <- cbind(compare, ratio = compare[, 2]/compare[, 
        3])
    colnames(compare) <- c("age", "sim_ratio", "true_ratio", 
        "ratio")
    the.factor <- min(compare[!is.na(compare[, "ratio"]), "ratio"])
    age.prob <- compare[(compare[!is.na(compare[, "ratio"]), 
        "ratio"] == min(compare[!is.na(compare[, "ratio"]), "ratio"])), 
        "age"]
    new.pedpop <- cbind(compare[, "age"], compare[, "true_ratio"], 
        rep(0, length = nrow(pedpop)))
    colnames(new.pedpop) <- c("age", "true_ratio", "counts")
    new.pedpop[new.pedpop[, "age"] == this.age, "counts"] <- the.factor * 
        pedpop[pedpop[, "age"] == this.age, "counts"]
    new.pedpop[, "counts"] <- new.pedpop[, "true_ratio"] * new.pedpop[new.pedpop[, 
        "age"] == this.age, "counts"]
    probs <- cbind(new.pedpop[, "age"], p = new.pedpop[, "counts"]/pedpop[, 
        "counts"])
    colnames(probs) <- c("age", "p")
    probs[probs[, "age"] == this.age, "p"] <- new.pedpop[probs[, 
        "age"] == this.age, "counts"]/old.value.this
    ped <- cbind(ped, age = (baseyr - ped[, "yob"]))
    ped <- cbind(ped, prob = probs[, "p"][match(ped[, "age"], 
        probs[, "age"])])
    ped <- ped[!is.na(ped[, "prob"]), ]
    ped[ped[, "prob"] > 1, "prob"] <- 1
    ped <- cbind(ped, p = rbinom(nrow(ped), 1, ped[, "prob"]))
    ped <- ped[ped[, "p"] == 1, ]
    ped <- ped[, c(1:7)]
    return(as.matrix(ped))
}
"draw.profile" <-
function (pop, year, realmale, realfem) 
{
    profile.male <- age.profile(pop, sex = 1, year)
    profile.fem <- age.profile(pop, sex = 2, year)
    profile.male <- profile.male[profile.male[, "age"] %in% realmale[, 
        "age"], ]
    realmale <- realmale[realmale[, "age"] %in% profile.male[, 
        "age"], ]
    profile.fem <- profile.fem[profile.fem[, "age"] %in% realfem[, 
        "age"], ]
    realfem <- realfem[realfem[, "age"] %in% profile.fem[, "age"], 
        ]
    realmale[, "counts"] <- (realmale[, "counts"]/sum(realmale[, 
        "counts"])) * 100
    realfem[, "counts"] <- (realfem[, "counts"]/sum(realfem[, 
        "counts"])) * 100
    profile.male[, "counts"] <- (profile.male[, "counts"]/sum(profile.male[, 
        "counts"])) * 100
    profile.fem[, "counts"] <- (profile.fem[, "counts"]/sum(profile.fem[, 
        "counts"])) * 100
    split.screen(c(1, 2))
    screen(1)
    par(cex = 0.7)
    plot(range(c(realmale[, "age"], profile.male[, "age"])), 
        range(c(realmale[, "counts"], profile.male[, "counts"])), 
        type = "n", xlab = "age", ylab = "percentage")
    title(paste("Age profiles for the real and simulated \n MALE population of year", 
        year))
    lines(realmale[, "age"], realmale[, "counts"], col = "red", 
        lwd = 3, lty = 3)
    lines(profile.male[, "age"], profile.male[, "counts"], col = "red", 
        lwd = 2)
    leg.txt <- c("Real population", "Simulated population")
    legend(0, 0.22, leg.txt, lwd = c(3, 2), lty = c(3, 1), col = c("red", 
        "red"), cex = 0.9, bty = "n")
    screen(2)
    par(cex = 0.7)
    plot(range(c(realfem[, "age"], profile.fem[, "age"])), range(c(realfem[, 
        "counts"], profile.fem[, "counts"])), type = "n", xlab = "age", 
        ylab = "percentage")
    title(paste("Age profiles for the real and simulated \n FEMALE population of year", 
        year))
    lines(realfem[, "age"], realfem[, "counts"], col = "darkblue", 
        lwd = 3, lty = 3)
    lines(profile.fem[, "age"], profile.fem[, "counts"], col = "darkblue", 
        lwd = 2)
    leg.txt <- c("Real population", "Simulated population")
    legend(0, 0.22, leg.txt, lwd = c(3, 2), lty = c(3, 1), col = c("darkblue", 
        "darkblue"), cex = 0.9, bty = "n")
}
"givebirth" <-
function (ped, childbc, nr, yr) 
{
    probfar <- ped[ped[, "yod"] == 0 & ped[, "sex"] == 1 & ped[, 
        "ID"] %in% setdiff(ped[ped[, "sex"] == 1, "ID"], unique(ped[ped[, 
        "f"] != 0, "f"])), , drop = FALSE]
    probfar <- cbind(probfar, age = (yr - probfar[, "yob"]))
    probmor <- ped[ped[, "yod"] == 0 & ped[, "yoi"] == 0 & ped[, 
        "sex"] == 2, , drop = FALSE]
    probmor <- cbind(probmor, age = (yr - probmor[, "yob"]))
    probmor <- cbind(probmor, rate = childbc[, "rate"][match(probmor[, 
        "age"], childbc[, "age"])])
    probmor <- probmor[!is.na(probmor[, "rate"]), ]
    if (nrow(probmor) == 0) 
        stop(paste("No potential mothers in the population. Calendar year", 
            yr))
    if (nrow(probmor[probmor[, "rate"] > 1, , drop = FALSE]) > 
        0) {
        warning("One (or more) probability of reproduction greater than 1 \n")
        probmor[probmor[, "rate"] > 1, "rate"] <- 1
    }
    already.mor <- probmor[probmor[, "ID"] %in% unique(ped[ped[, 
        "m"] != 0, "m"]), c("ID", "rate")]
    already.couple <- cbind(m = already.mor[, "ID"], rate = already.mor[, 
        "rate"], f = ped[, "f"][match(already.mor[, "ID"], ped[, 
        "m"])])
    already.couple <- already.couple[!is.na(already.couple[, 
        "f"]), ]
    already.couple <- cbind(already.couple, yod.f = ped[, "yod"][match(already.couple[, 
        "f"], ped[, "ID"])])
    already.couple <- already.couple[!is.na(already.couple[, 
        "yod.f"]) & already.couple[, "yod.f"] == 0, ]
    already.couple <- cbind(already.couple, chances = rbinom(nrow(already.couple), 
        1, already.couple[, "rate"]))
    first.timer <- probmor[probmor[, "ID"] %in% setdiff(probmor[, 
        "ID"], already.mor[, "ID"]), ]
    first.timer <- cbind(m = first.timer[, "ID"], age = first.timer[, 
        "age"], chances = rbinom(nrow(first.timer), 1, first.timer[, 
        "rate"]), f = rep(0, length = nrow(first.timer)))
    first.timer <- first.timer[first.timer[, "chances"] == 1, 
        ]
    ages <- unique(first.timer[, "age"])
    for (a in 1:length(ages)) {
        if (nrow(first.timer[first.timer[, "age"] == ages[a], 
            , drop = FALSE]) <= nrow(probfar[probfar[, "age"] %in% 
            c((ages[a] - 1):(ages[a] + 4)), , drop = FALSE])) {
            first.timer[first.timer[, "age"] == ages[a], "f"] <- sample(probfar[probfar[, 
                "age"] %in% c((ages[a] - 1):(ages[a] + 4)), "ID", 
                drop = FALSE], size = nrow(first.timer[first.timer[, 
                "age"] == ages[a], , drop = FALSE]))
        }
        else {
            first.timer[first.timer[, "m"] %in% sample(first.timer[first.timer[, 
                "age"] == ages[a] & first.timer[, "f"] == 0, 
                "m"], size = nrow(probfar[probfar[, "age"] %in% 
                c((ages[a] - 1):(ages[a] + 4)), , drop = FALSE])), 
                "f"] <- probfar[probfar[, "age"] %in% c((ages[a] - 
                1):(ages[a] + 4)), "ID"]
        }
        probfar <- probfar[probfar[, "ID"] %in% setdiff(probfar[, 
            "ID"], first.timer[, "f"]), ]
    }
    mor.far <- rbind(already.couple[already.couple[, "chances"] == 
        1, c("m", "f"), drop = FALSE], first.timer[first.timer[, 
        "f"] != 0, c("m", "f"), drop = FALSE])
    mor.far <- cbind(mor.far, sex = rbinom(nrow(mor.far), 1, 
        0.5))
    ped <- rbind(ped, cbind(c((nr + 1):(nr + nrow(mor.far))), 
        rep(yr, length = nrow(mor.far)), ifelse(mor.far[, "sex"] == 
            1, 2, 1), mor.far[, c("m", "f")], rep(0, length = nrow(mor.far)), 
        rep(0, length = nrow(mor.far))))
    return(ped)
}
"plot.poplab" <-
function (x, option, population.fem, population.male, year, folder, 
    ...) 
{
    if (missing(option)) {
        option <- "current"
    }
    if (option == "current") 
        pop.t <- x[[2]]
    if (option == "base") {
        if (is.null(dim(x)[2])) 
            pop.t <- x[[1]]
        if (!is.null(dim(x)[2])) 
            pop.t <- x
    }
    realfem <- readin(f.name = population.fem, folder)
    realmale <- readin(f.name = population.male, folder)
    if ((!(year %in% as.integer(colnames(realfem)[-1]))) | (!(year %in% 
        as.integer(colnames(realmale)[-1])))) 
        stop(paste("No female or male population counts for the year", 
            year, "in data files"))
    realmale <- realmale[, c(1, which(colnames(realmale) == year, 
        arr.ind = TRUE))]
    colnames(realmale) <- c("age", "counts")
    realfem <- realfem[, c(1, which(colnames(realfem) == year, 
        arr.ind = TRUE))]
    colnames(realfem) <- c("age", "counts")
    draw.profile(pop.t, year, realmale, realfem)
}
"print.poplab" <-
function (x, option, year, folder, ...) 
{
    require(MASS)
    if (option == "current") {
        pop.t <- x[[2]]
        colnames(pop.t) <- c("ID", "yob", "sex", "m", "f", "yod", 
            "yoi")
        pop.t <- pop.t[pop.t[, "yob"] <= year, ]
    }
    if (option == "base") {
        if (is.null(dim(x)[2])) 
            pop.t <- x[[1]]
        if (!is.null(dim(x)[2])) 
            pop.t <- x
    }
    where <- switch(option, base = paste("base_pop_", year, ".txt", 
        sep = ""), current = paste("simpop_endYr_", year, ".txt", 
        sep = ""))
    write.matrix(pop.t, file = file.path(folder, where), sep = " ")
    cat("Data saved to", folder, "folder as", where, "file \n")
}
"readin" <-
function (f.name, folder) 
{
    name <- scan(file = file.path(folder, f.name), quiet = TRUE, 
        what = "character", nlines = 1)
    pop <- matrix(scan(file = file.path(folder, f.name), skip = 1, 
        quiet = TRUE), ncol = length(name), byrow = TRUE)
    colnames(pop) <- name
    pop
}
"retrieve.parents" <-
function (pop, bigpop) 
{
    pop <- pop[, c("ID", "yob", "sex", "m", "f", "yod", "yoi")]
    parents <- unique(c(pop[, "m"], pop[, "f"]))
    parents <- parents[parents %in% setdiff(parents, pop[, "ID"])]
    notaround <- bigpop[match(parents, bigpop[, "ID"]), ]
    colnames(notaround) <- c("ID", "yob", "sex", "m", "f", "yod", 
        "yoi")
    notaround[is.na(notaround[, "ID"]), "ID"] <- 0
    notaround[is.na(notaround[, "yob"]), "yob"] <- 0
    notaround[is.na(notaround[, "sex"]), "sex"] <- 0
    notaround[is.na(notaround[, "m"]), "m"] <- 0
    notaround[is.na(notaround[, "f"]), "f"] <- 0
    notaround[is.na(notaround[, "yod"]), "yod"] <- 0
    notaround[is.na(notaround[, "yoi"]), "yoi"] <- 0
    notaround <- notaround[notaround[, "ID"] != 0 & notaround[, 
        "yob"] != 0, ]
    return(notaround)
}
"simped" <-
function (baseyear, basehealth, basefamrisk, sex.a, basetotal, 
    warmuptime, simyears, endyear, healthy, famrisk, mortratio, 
    d.mod, fam.rel, print.option, population.fem, population.male, 
    mortality.fem, mortality.male, fertility, incidence, seed, 
    folder, name.base) 
{
    if (missing(basehealth)) {
        basehealth <- 1
    }
    if (missing(basefamrisk)) {
        basefamrisk <- 1
    }
    if (missing(sex.a)) {
        sex.a <- 2
    }
    if (missing(basetotal)) {
        basetotal <- 500
    }
    if (missing(warmuptime)) {
        warmuptime <- 100
    }
    if (missing(endyear)) {
        endyear <- baseyear + simyears - 1
    }
    if (missing(baseyear)) {
        baseyear <- endyear - simyears + 1
    }
    if (missing(healthy)) {
        healthy <- 1
    }
    if (missing(famrisk)) {
        famrisk <- 1
    }
    if (missing(mortratio)) {
        mortratio <- 1
    }
    if (missing(print.option)) {
        print.option <- FALSE
    }
    if (missing(seed)) {
        seed <- NULL
    }
    if (missing(name.base)) {
        name.base <- ""
    }
    if ((!is.null(seed)) & (!is.character(seed))) 
        set.seed(seed)
    if (basetotal <= 0) 
        stop("Something wrong with parameter basetotal. See help(simped)")
    if (warmuptime <= 0) 
        stop("Something wrong with parameter warmuptime. See help(simped)")
    fert <- readin(f.name = fertility, folder)
    femmort <- readin(f.name = mortality.fem, folder)
    malemort <- readin(f.name = mortality.male, folder)
    warn.old <- options()$warn
    options(warn = -1)
    if (name.base == "") {
        name <- try(scan(file = file.path(folder, paste("base_pop_", baseyear, 
            ".txt", sep = "")), quiet = TRUE, what = "character", 
            nlines = 1), silent = TRUE)
    }
    else {
        name <- try(scan(file = file.path(folder, paste("base_pop_", baseyear, 
            "_", name.base, ".txt", sep = "")), quiet = TRUE, 
            what = "character", nlines = 1), silent = TRUE)
    }
    if (!(name[1] %in% c("ID", "yob", "sex", "m", "f", "yod", 
        "yoi"))) {
        if (basehealth == 0) {
            if (length(basefamrisk) > 1) {
                basefamrisk <- matrix(data = basefamrisk, byrow = FALSE, 
                  ncol = 2)
                colnames(basefamrisk) <- c("age", "val")
                if (!(d.mod %in% c("agesprr", "agespor"))) 
                  stop(paste("Invalid risk value specification for the chosen disease model"))
            }
            else {
                if (!(d.mod %in% c("rr", "or"))) 
                  stop(paste("Invalid risk value specification for the chosen disease model"))
            }
            if (!(d.mod %in% c("agespor", "agesprr", "or", "rr"))) 
                stop(paste("Invalid disease model specification for the baseline population"))
        }
        basepop <- create.baseline.complete(baseyear, healthy = basehealth, 
            risk = basefamrisk, sex.a = sex.a, mortratio = mortratio, 
            base.scale = basetotal, runintime = warmuptime, d.mod = d.mod, 
            fam.rel = fam.rel, print.option = print.option, population.fem = population.fem, 
            population.male = population.male, mortality.fem = mortality.fem, 
            mortality.male = mortality.male, fertility = fertility, 
            incidence = incidence, seed = seed, folder = folder)
        plot.poplab(basepop, "base", population.fem, population.male, 
            baseyear, folder)
    }
    if (name[1] %in% c("ID", "yob", "sex", "m", "f", "yod", "yoi")) {
        if (name.base == "") {
            basepop <- readin(f.name = paste("base_pop_", baseyear, 
                ".txt", sep = ""), folder)
        }
        else {
            basepop <- readin(f.name = paste("base_pop_", baseyear, 
                "_", name.base, ".txt", sep = ""), folder)
        }
        cat("Baseline population was read from", folder, "folder \n")
        class(basepop) <- "poplab"
        plot.poplab(basepop, "base", population.fem, population.male, 
            baseyear, folder)
    }
    options(warn = warn.old)
    pedigree <- basepop
    todel <- which(pedigree[, "yob"] == baseyear)
    if (length(todel) != 0) 
        pedigree <- pedigree[-todel, ]
    counter <- max(pedigree[, "ID"])
    cat("Start simulating pedigree from calendar year", baseyear, 
        "to calendar year", endyear, ": \n")
    cat("Baseline population:", nrow(pedigree), "individuals \n")
    if (endyear > max(as.numeric(colnames(fert)[-1]))) 
        warning(paste("Fertility rates of the year", max(as.numeric(colnames(fert)[-1])), 
            "will be used for the calendar years", max(as.numeric(colnames(fert)[-1])) + 
                1, "-", endyear), call. = FALSE)
    if (endyear > max(as.numeric(colnames(femmort)[-1]))) 
        warning(paste("Female mortality rates of the year", max(as.numeric(colnames(femmort)[-1])), 
            "will be used for the calendar years", max(as.numeric(colnames(femmort)[-1])) + 
                1, "-", endyear), call. = FALSE)
    if (endyear > max(as.numeric(colnames(malemort)[-1]))) 
        warning(paste("Male mortality rates of the year", max(as.numeric(colnames(malemort)[-1])), 
            "will be used for the calendar years", max(as.numeric(colnames(malemort)[-1])) + 
                1, "-", endyear), call. = FALSE)
    if (healthy == 0) {
        inc <- readin(f.name = incidence, folder)
        if (endyear > max(as.numeric(colnames(inc)[-1]))) 
            warning(paste("Incidence rates of the year", max(as.numeric(colnames(inc)[-1])), 
                "will be used for the calendar years", max(as.numeric(colnames(inc)[-1])) + 
                  1, "-", endyear), call. = FALSE)
        if (length(famrisk) > 1) {
            famrisk <- matrix(data = famrisk, byrow = FALSE, 
                ncol = 2)
            colnames(famrisk) <- c("age", "val")
            if (!(d.mod %in% c("agesprr", "agespor"))) 
                stop(paste("Invalid risk value specification for the chosen disease model"))
        }
        else {
            if (!(d.mod %in% c("rr", "or"))) 
                stop(paste("Invalid risk value specification for the chosen disease model"))
        }
        if (!(d.mod %in% c("agespor", "agesprr", "or", "rr"))) 
            stop(paste("Invalid disease model specification for the evolved population"))
    }
    z <- baseyear
    while (z <= endyear) {
        ind <- ifelse(z > max(as.numeric(colnames(fert)[-1])), 
            1 + which.max(as.numeric(colnames(fert)[-1])), which(colnames(fert) == 
                z, arr.ind = TRUE))
        childbc <- fert[, c(1, ind)]
        childbc <- childbc[childbc[, 2] != 0, ]
        colnames(childbc) <- c("age", "rate")
        indmort.f <- ifelse(z > max(as.numeric(colnames(femmort)[-1])), 
            1 + which.max(as.numeric(colnames(femmort)[-1])), 
            which(colnames(femmort) == z, arr.ind = TRUE))
        indmort.m <- ifelse(z > max(as.numeric(colnames(malemort)[-1])), 
            1 + which.max(as.numeric(colnames(malemort)[-1])), 
            which(colnames(malemort) == z, arr.ind = TRUE))
        if (nrow(femmort) != nrow(malemort)) 
            stop("Different age structure in the male and female mortality files. Organize files and try again.")
        tempmort <- cbind(femmort[, c(1, indmort.f)], malemort[, 
            indmort.m])
        colnames(tempmort) <- c("age", "femrate", "malerate")
        pedigree <- givebirth(ped = pedigree, childbc = childbc, 
            nr = counter, yr = z)
        counter <- max(pedigree[, "ID"])
        pedigree <- assigndeath(ped = pedigree, tempmort = tempmort, 
            risk = mortratio, yr = z)
        if (healthy == 0) {
            ind.i <- ifelse(z > max(as.numeric(colnames(inc)[-1])), 
                1 + which.max(as.numeric(colnames(inc)[-1])), 
                which(colnames(inc) == z, arr.ind = TRUE))
            tempinc <- inc[, c(1, ind.i)]
            tempinc <- tempinc[tempinc[, 2] != 0, , drop = FALSE]
            colnames(tempinc) <- c("age", "rate")
            pedigree <- assigncancer(ped = pedigree, inc = tempinc, 
                risk = famrisk, yr = z, sex = sex.a, type = d.mod, 
                fam.rel = fam.rel)
        }
        z <- z + 1
    }
    cat("Finished simulating the population for calendar time", 
        baseyear, "-", endyear, "\n")
    class(pedigree) <- "poplab"
    res <- list(basepop, pedigree)
    class(res) <- "poplab"
    if (print.option) 
        print.poplab(res, "current", endyear, folder)
    return(res)
}
