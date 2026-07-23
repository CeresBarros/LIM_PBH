# syntax=docker/dockerfile:1

######################
## install git and get files before creating the image

FROM ubuntu:latest AS builder

## install git and get git repo with installation scripts

RUN apt-get update && \
apt-get upgrade -y && \
apt-get install -y git

WORKDIR /tmp

RUN git clone https://github.com/rocker-org/rocker-versioned2.git rocker-versioned2/

## edit a few of the scripts

## edit *.sh to patch R config file to change required libcurl version
## https://unix.stackexchange.com/questions/758026/error-libcurl-7-28-0-library-and-headers-are-required-with-support-for-https
RUN <<EOF
awk '
  { print }
  /^cd R-\*\// {
    print "\nif [ \"$R_VERSION\" == \"4.1.3\" ]; then"
    print "  awk \047/> 7/ { c = 1 } !/> 7/ && c { print(\"  exit(0);\"); c = 0; next; } 1\047 configure > configure.new && mv configure.new configure"
    print "fi"
    print "\nchmod +x configure"
  }' rocker-versioned2/scripts/install_R_source.sh > temp.sh
mv temp.sh rocker-versioned2/scripts/install_R_source.sh
chmod +x rocker-versioned2/scripts/install_R_source.sh
EOF


## edit wget when downloading various software pieces
RUN awk '{gsub(/wget\ \"https/, "wget --no-check-certificate \"http")}1' rocker-versioned2/scripts/install_R_source.sh > temp.sh && \
mv temp.sh rocker-versioned2/scripts/install_R_source.sh && \
chmod +x rocker-versioned2/scripts/install_R_source.sh

RUN <<EOF
awk '{gsub(/wget -P/, "wget --no-check-certificate -P")}1' rocker-versioned2/scripts/install_s6init.sh > temp.sh
mv temp.sh rocker-versioned2/scripts/install_s6init.sh
awk '{gsub(/https/, "http")}1' rocker-versioned2/scripts/install_s6init.sh > temp.sh
mv temp.sh rocker-versioned2/scripts/install_s6init.sh
chmod +x rocker-versioned2/scripts/install_s6init.sh
EOF

RUN awk '{gsub(/wget\ \"https/, "wget --no-check-certificate \"http")}1' rocker-versioned2/scripts/install_rstudio.sh > temp.sh && \
mv temp.sh rocker-versioned2/scripts/install_rstudio.sh && \
chmod +x rocker-versioned2/scripts/install_rstudio.sh

RUN awk '{gsub(/wget\ \"\$PANDOC_DL_URL/, "wget --no-check-certificate \"$PANDOC_DL_URL")}1' rocker-versioned2/scripts/install_pandoc.sh > temp.sh && \
mv temp.sh rocker-versioned2/scripts/install_pandoc.sh && \
chmod +x rocker-versioned2/scripts/install_pandoc.sh

RUN awk '{gsub(/wget\ \"\$QUARTO_DL_URL/, "wget --no-check-certificate \"$QUARTO_DL_URL")}1' rocker-versioned2/scripts/install_quarto.sh > temp.sh && \
mv temp.sh rocker-versioned2/scripts/install_quarto.sh && \
chmod +x rocker-versioned2/scripts/install_quarto.sh

######################
## Now create image

FROM docker.io/library/ubuntu:latest

ENV R_VERSION="4.1.3"
ENV R_HOME="/usr/local/lib/R"
ENV TZ="Etc/UTC"
ENV HOME="/home/rstudio"
ENV USER="rstudio"

## install necessary software
RUN apt-get update && apt-get install -y openssh-client \
  git \
  vim \
  less \
  nano \
  htop \
  rsync \
  libglpk40 \
  libgdal-dev \
  gdal-bin \
  proj-bin \
  libproj-dev \
  libpng-dev \
  libudunits2-dev \
  libgeos-dev \
  libsqlite3-dev \
  locate \
  r-cran-littler \
  libabsl-dev \
  libmagick++-dev

## user configs
RUN <<EOF
mkdir $HOME
groupadd -g 1003 "$USER"
useradd -m -u 1003 -g "$USER" "$USER"
useradd sudo "$USER"
## change ownership of homedir
chown -R "$USER":"$USER" "$HOME"
usermod -aG sudo "$USER"
EOF

## need to configure a few things.
RUN <<EOF
export LD_LIBRARY_PATH="~/bin/gdal/lib:$LD_LIBRARY_PATH"
export LD_LIBRARY_PATH="/usr/local/lib:$LD_LIBRARY_PATH"
ldconfig -c "echo '/usr/local/lib' >> /etc/ld.so.conf.d/R-dependencies-x86_64.conf"
ln -s libgdal.so /usr/lib/x86_64-linux-gnu/libgdal.so.26
ln -s libproj.so /usr/lib/x86_64-linux-gnu/libproj.so.15
ln -s /usr/lib/x86_64-linux-gnu/libproj.so /usr/lib/libproj.so
ln -s /usr/lib/x86_64-linux-gnu/libproj.so /usr/lib/libproj.so.15
ln -s /usr/lib/x86_64-linux-gnu/libMagick++-6.Q16.so /usr/lib/libMagickCore-6.Q16.so.6
mkdir -p /usr/local/gdal/share
ln -s /usr/share/gdal /usr/local/gdal/share/gdal
export GDAL_DATA="/usr/local/gdal/share/gdal"
mkdir -p /usr/local/gdal/bin/
ln -s /usr/bin/gdal* /usr/local/gdal/bin/
EOF

COPY --from=builder /tmp/rocker-versioned2/scripts/install_R_source.sh /rocker_scripts/install_R_source.sh
RUN /rocker_scripts/install_R_source.sh

ENV CRAN="http://cloud.r-project.org"
ENV LANG=en_US.UTF-8

COPY --from=builder /tmp/rocker-versioned2/scripts/bin/ /rocker_scripts/bin/
COPY --from=builder /tmp/rocker-versioned2/scripts/setup_R.sh /rocker_scripts/setup_R.sh

RUN <<EOF
if grep -q "1000" /etc/passwd; then
userdel --remove "$(id -un 1000)";
fi
/rocker_scripts/setup_R.sh
EOF

ENV S6_VERSION="v2.1.0.2"
ENV RSTUDIO_VERSION="2025.05.1+513"
ENV DEFAULT_USER="rstudio"

COPY --from=builder /tmp/rocker-versioned2/scripts/install_rstudio.sh /rocker_scripts/install_rstudio.sh
COPY --from=builder /tmp/rocker-versioned2/scripts/install_s6init.sh /rocker_scripts/install_s6init.sh
COPY --from=builder /tmp/rocker-versioned2/scripts/default_user.sh /rocker_scripts/default_user.sh
COPY --from=builder /tmp/rocker-versioned2/scripts/init_set_env.sh /rocker_scripts/init_set_env.sh
COPY --from=builder /tmp/rocker-versioned2/scripts/init_userconf.sh /rocker_scripts/init_userconf.sh
COPY --from=builder /tmp/rocker-versioned2/scripts/pam-helper.sh /rocker_scripts/pam-helper.sh
RUN /rocker_scripts/install_rstudio.sh

EXPOSE 8787
CMD ["/init"]

COPY --from=builder /tmp/rocker-versioned2/scripts/install_pandoc.sh /rocker_scripts/install_pandoc.sh
RUN /rocker_scripts/install_pandoc.sh

COPY --from=builder /tmp/rocker-versioned2/scripts/install_quarto.sh /rocker_scripts/install_quarto.sh
RUN /rocker_scripts/install_quarto.sh

COPY --from=builder /tmp/rocker-versioned2/scripts /rocker_scripts
